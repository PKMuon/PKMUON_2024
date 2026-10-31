"""Read-only extraction of Altium compound streams. Never executes PCB content."""
import struct,json,hashlib,re
from pathlib import Path
import sys
def streams(path):
 b=Path(path).read_bytes();assert b[:8]==bytes.fromhex('d0cf11e0a1b11ae1')
 u=lambda off:struct.unpack_from('<I',b,off)[0]
 size=1<<struct.unpack_from('<H',b,30)[0];mini=1<<struct.unpack_from('<H',b,32)[0]
 sec=lambda n:b[(n+1)*size:(n+2)*size]
 dif=list(struct.unpack_from('<109I',b,76));n=u(68)
 for _ in range(u(72)):
  a=struct.unpack('<'+'I'*(size//4),sec(n));dif.extend(a[:-1]);n=a[-1]
 fat=[]
 for n in dif:
  if n<0xfffffffa:fat.extend(struct.unpack('<'+'I'*(size//4),sec(n)))
 def chain(n,table,get):
  out=[];seen=set()
  while n<0xfffffffa:
   assert n not in seen;seen.add(n);out.append(get(n));n=table[n]
  return b''.join(out)
 raw=chain(u(48),fat,sec);entries=[]
 for off in range(0,len(raw),128):
  a=raw[off:off+128];ln=struct.unpack_from('<H',a,64)[0]
  entries.append(dict(name=a[:max(0,ln-2)].decode('utf-16le',errors='replace'),typ=a[66],left=struct.unpack_from('<I',a,68)[0],right=struct.unpack_from('<I',a,72)[0],child=struct.unpack_from('<I',a,76)[0],start=struct.unpack_from('<I',a,116)[0],length=struct.unpack_from('<Q',a,120)[0]))
 root=entries[0];small=chain(root['start'],fat,sec)[:root['length']]
 mfraw=chain(u(60),fat,sec) if u(64) else b'';mf=list(struct.unpack('<'+'I'*(len(mfraw)//4),mfraw))
 result={}
 def walk(i,parent):
  if i>=len(entries):return
  e=entries[i];walk(e['left'],parent);name=parent+[e['name']]
  if e['typ']==1:walk(e['child'],name)
  if e['typ']==2:
   data=chain(e['start'],mf,lambda n:small[n*mini:(n+1)*mini]) if e['length']<u(56) else chain(e['start'],fat,sec)
   result['/'.join(name)]=data[:e['length']]
  walk(e['right'],parent)
 walk(root['child'],[]);return result
def extract(path,out):
 out=Path(out);out.mkdir(parents=True,exist_ok=True);ss=streams(path);manifest={}
 for name,data in ss.items():
  dest=out/(name.replace('/','__')+'.bin');dest.write_bytes(data)
  strings=[s.decode('gb18030',errors='replace') for s in re.findall(rb'[\x20-\x7e\x80-\xff]{12,}',data)]
  manifest[name]={'bytes':len(data),'strings':strings}
 (out/'streams.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
 print('\n'.join(f'{k}: {len(v)}' for k,v in ss.items()))
 return ss
if __name__=='__main__':extract(sys.argv[1],sys.argv[2])
