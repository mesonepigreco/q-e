#!/usr/bin/env python3
"""Add explicit frozen-core tau to a matching open PAW UPF from PAW-XML.

The XML radial core fields are Y00 coefficients. UPF fields are physical
scalar densities. This converter does not alter projectors, waves or potentials.
Requires numpy/scipy. Both input files must describe the SAME atomic dataset.
"""
import argparse,hashlib,json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import simpson

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('xml',type=Path);p.add_argument('upf',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
 x=ET.parse(a.xml).getroot();u=ET.parse(a.upf).getroot();atom=x.find('atom');h=u.find('PP_HEADER')
 if h.attrib['pseudo_type']!='PAW' or atom.attrib['symbol']!=h.attrib['element']:p.error('Species/PAW mismatch')
 if abs(float(atom.attrib['valence'])-float(h.attrib['z_valence']))>1.e-8:p.error('Valence mismatch')
 r=np.fromstring(x.find('radial_grid/values').text,sep=' ');q=np.fromstring(u.find('PP_MESH/PP_R').text,sep=' ')
 density=np.fromstring(x.find('ae_core_density').text,sep=' ')/np.sqrt(4*np.pi)
 charge=float(4*np.pi*simpson(density*r*r,x=r));expected=float(atom.attrib['core'])
 if abs(charge-expected)>1.e-5:p.error('XML core normalization failed')
 upf_density=np.fromstring(u.find('PP_PAW/PP_AE_NLCC').text,sep=' ')
 sample=CubicSpline(r,density)(q)
 relative=float(np.max(np.abs(sample-upf_density))/max(np.max(np.abs(upf_density)),1.e-30))
 if relative>.005:p.error('AE core densities are inconsistent; dataset pairing must be checked')
 # Compare both AE and PS partial waves in addition to the core density.
 # PAW-XML stores R(r); UPF stores u(r)=r*R(r). Never mix tau from a
 # merely similar atom or functional with a different PAW operator.
 wave_errors=[]
 states=x.findall('valence_states/state')
 for family,tag in [('ae_partial_wave','PP_AEWFC.'),('pseudo_partial_wave','PP_PSWFC.')]:
  waves=x.findall(family)
  if len(waves)!=int(h.attrib['number_of_proj']):p.error('Projector count mismatch')
  for i,(state,wave) in enumerate(zip(states,waves),1):
   node=u.find('PP_FULL_WFC/'+tag+str(i))
   if node is None or int(node.attrib['l'])!=int(state.attrib['l']):p.error('Angular channel mismatch')
   values=np.fromstring(wave.text,sep=' ');target=np.fromstring(node.text,sep=' ')
   err=float(np.max(np.abs(CubicSpline(r,values)(q)*q-target))/max(np.max(np.abs(target)),1.e-30))
   wave_errors.append(err)
   if err>1.e-5:p.error('Partial waves do not describe the same atomic dataset')
 paw=u.find('PP_PAW');paw.set('tb09_core_tau','T')
 for source,target in [('ae_core_kinetic_energy_density','PP_AE_TAUCORE'),('pseudo_core_kinetic_energy_density','PP_PS_TAUCORE')]:
  node=x.find(source)
  if node is None:p.error('Missing XML field '+source)
  values=np.fromstring(node.text,sep=' ')/np.sqrt(4*np.pi)
  if np.any(~np.isfinite(values)) or np.min(values)<-1.e-12:p.error('Invalid positive kinetic density')
  result=CubicSpline(r,values)(q)
  if result.min()<-1.e-8*max(1,result.max()):p.error('Interpolation introduced negative tau')
  result=np.maximum(result,0)
  old=paw.find(target)
  if old is not None:paw.remove(old)
  v=ET.SubElement(paw,target,{'size':str(len(q)),'units':'Ha/bohr^3; tau=0.5*sum|grad(phi)|^2'})
  v.tail='\n'
  v.text='\n'+'\n'.join(' '.join(f'{z:.17e}' for z in result[j:j+3]) for j in range(0,len(result),3))+'\n'
 ET.ElementTree(u).write(a.output,encoding='unicode')
 manifest={'xml_sha256':hashlib.sha256(a.xml.read_bytes()).hexdigest(),'upf_sha256':hashlib.sha256(a.upf.read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(a.output.read_bytes()).hexdigest(),'AE_core_electrons':charge,'core_density_pair_relative_error':relative,'partial_wave_pair_max_relative_error':max(wave_errors),'tau_convention':'physical scalar Hartree density; XML Y00 divided by sqrt(4*pi)'}
 a.output.with_suffix(a.output.suffix+'.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest,indent=2))
if __name__=='__main__':main()
