"""Decode supported Altium v6 records with strict framing/count checks."""
from pathlib import Path
import struct,json,hashlib,csv
from inspect_pcb import streams
SOURCE=Path('C:/Users/Brohuy/Desktop/杂类/新2维读出板-悬浮延迟块引脚改接地(1).PCB')
OUT=Path(__file__).resolve().parents[1]/'results/strip_readout_20261030';OUT.mkdir(parents=True,exist_ok=True)
s=streams(SOURCE);scale=.0254/10000
def props(b):return dict(x.split('=',1) for x in b.rstrip(b'\0').decode('gb18030',errors='replace').split('|') if '=' in x)
def blocks(b):
 o=0;out=[]
 while o<len(b):
  n=struct.unpack_from('<I',b,o)[0];o+=4;assert o+n<=len(b);out.append(b[o:o+n]);o+=n
 return out
def count(k):return struct.unpack('<I',s[k+'/Header'])[0]
nets=[props(b) for b in blocks(s['Nets6/Data'])];components=[props(b) for b in blocks(s['Components6/Data'])]
assert len(nets)==count('Nets6') and len(components)==count('Components6')
def common(b):
 net=struct.unpack_from('<H',b,3)[0];comp=struct.unpack_from('<H',b,7)[0]
 return dict(layer=b[0],net_id=net,net=nets[net]['NAME'] if net!=65535 else None,component_id=comp)
def records(k,typ):
 b=s[k+'/Data'];out=[];o=0
 while o<len(b):
  assert b[o]==typ;o+=1;n=struct.unpack_from('<I',b,o)[0];o+=4;out.append(b[o:o+n]);o+=n
 assert o==len(b) and len(out)==count(k);return out
fills=[]
for b in records('Fills6',6):
 d=common(b);d.update(zip(['x1_mm','y1_mm','x2_mm','y2_mm'],[a*scale for a in struct.unpack_from('<4i',b,13)]));d['rotation_deg']=struct.unpack_from('<d',b,29)[0];fills.append(d)
tracks=[]
for b in records('Tracks6',4):
 d=common(b);d.update(zip(['x1_mm','y1_mm','x2_mm','y2_mm','width_mm'],[a*scale for a in struct.unpack_from('<5i',b,13)]));tracks.append(d)
pads=[];b=s['Pads6/Data'];o=0
while o<len(b):
 assert b[o]==2;o+=1;bb=[]
 for _ in range(6):
  n=struct.unpack_from('<I',b,o)[0];o+=4;bb.append(b[o:o+n]);o+=n
 d=common(bb[4]);d['pin']=bb[0][1:].decode('ascii');d.update(zip(['x_mm','y_mm','size_x_mm','size_y_mm'],[a*scale for a in struct.unpack_from('<4i',bb[4],13)]));pads.append(d)
assert len(pads)==count('Pads6')
board=props(blocks(s['Board6/Data'])[0]);stack={k:v for k,v in board.items() if k in ['LAYER1COPTHICK','LAYER1DIELCONST','LAYER1DIELHEIGHT','LAYER1DIELMATERIAL','TOPCONST','TOPHEIGHT','BOTTOMCONST','BOTTOMHEIGHT']}
data=dict(source=str(SOURCE),sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),units='mm',status='decoded copper primitives; not yet complete FEM geometry; arcs/vias/regions/polygons retained in raw streams',stackup_declared=stack,nets=[n['NAME'] for n in nets],components=components,fills=fills,tracks=tracks,pads=pads)
(Path(__file__).resolve().parents[2]/'config/pcb_geometry.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
with (OUT/'pad_nets.csv').open('w',newline='',encoding='utf-8-sig') as f:
 w=csv.DictWriter(f,fieldnames=list(pads[0]));w.writeheader();w.writerows(pads)
print(json.dumps({'counts':{k:len(data[k]) for k in ['nets','components','fills','tracks','pads']},'stackup':stack,'net_names':data['nets']},ensure_ascii=False))
