#Derived from Sunfish by Thomas Ahle (MIT)
import sys
V={'P':100,'N':320,'B':330,'R':500,'Q':900,'K':20000}
D={'N':[(-2,-1),(-2,1),(-1,-2),(-1,2),(1,-2),(1,2),(2,-1),(2,1)],'B':[(-1,-1),(-1,1),(1,-1),(1,1)],'R':[(-1,0),(1,0),(0,-1),(0,1)]}
D['Q']=D['K']=D['B']+D['R']
def G(b,tn,ep):
 w=tn=='w';ms=[]
 for i,p in enumerate(b):
  if p=='.'or p.isupper()!=w:continue
  r,f=i//8,i%8;u=p.upper()
  if u=='P':
   d=-8 if w else 8;t=i+d;a=lambda t:[(i,t,x)for x in('QRBN'if w else'qrbn')]if t//8 in(0,7)else[(i,t,'')]
   if b[t]=='.':
    ms+=a(t)
    if r==(6 if w else 1)and b[i+2*d]=='.':ms.append((i,i+2*d,''))
   for df in(-1,1):
    t=i+d+df
    if 0<=f+df<8 and(b[t]!='.'and b[t].isupper()!=w or t==ep):ms+=a(t)
  else:
   for dr,df in D[u]:
    x,y=r+dr,f+df
    while 0<=x<8 and 0<=y<8:
     t=x*8+y
     if b[t]=='.':ms.append((i,t,''))
     else:
      if b[t].isupper()!=w:ms.append((i,t,''))
      break
     if u in'NK':break
     x+=dr;y+=df
 return ms
def M(b,tn,cas,ep,m):
 f,t,pr=m;b=b[:];p=b[f];w=p.isupper();b[f]='.'
 if p.upper()=='P'and t==ep and b[t]=='.':b[t+(8 if w else-8)]='.'
 b[t]=pr or p;n=cas
 if p=='K':n=n.replace('K','').replace('Q','')
 if p=='k':n=n.replace('k','').replace('q','')
 for s,c in zip((63,56,7,0),'KQkq'):
  if f==s or t==s:n=n.replace(c,'')
 if p.upper()=='K'and abs(t-f)==2:
  rf,rt=(t+1,t-1)if t>f else(t-2,t+1);b[rt]=b[rf];b[rf]='.'
 return b,('b'if w else'w'),n,(t if p.upper()=='P'and abs(t-f)==16 else None)
def A(b,w):
 k=b.index('K'if w else'k');return any(m[1]==k for m in G(b,'b'if w else'w',None))
def L(b,tn,cas,ep):
 w=tn=='w';ms=[m for m in G(b,tn,ep)if not A(M(b,tn,cas,ep,m)[0],w)]
 i=60 if w else 4;K='K'if w else'k'
 if b[i]==K and not A(b,w):
  for sg,rr,rk in((1,i+3,K),(-1,i-4,chr(ord(K)+6))):
   if rk in cas and b[rr]==chr(ord(K)+7)and all(b[i+d]=='.'for d in range(sg,(3+(sg<0))*sg,sg)):
    nb=b[:];nb[i]='.';nb[i+sg]=K
    if not A(nb,w):
     nb[i+sg]='.';nb[i+2*sg]=K
     if not A(nb,w):ms.append((i,i+2*sg,''))
 ms.sort(key=lambda m:V.get(b[m[0]],0)-4096*V.get(b[m[1]].upper(),0));return ms
def E(b):
 s=0
 for i,p in enumerate(b):
  if p=='.':continue
  u=p.upper();r=i//8;v=V[u]+((6-r if p.isupper()else r-1)*7 if u=='P'else int(6-abs(i%8-3.5)-abs(r-3.5))*3 if u in'NBQ'else 0)
  s+=v if p.isupper()else-v
 return s
def S(b,tn,cas,ep,d,a,be,q=0):
 if q:
  st=(2*(tn=='w')-1)*E(b)
  if st>=be:return be
  a=max(a,st)
 ms=L(b,tn,cas,ep)
 if not ms:return-99999 if A(b,tn=='w')else 0
 if d==0 and not q:return S(b,tn,cas,ep,0,a,be,1)
 for m in ms:
  if q and b[m[1]]=='.'and m[1]!=ep and not m[2]:continue
  nb,nt,nc,ne=M(b,tn,cas,ep,m);sc=-S(nb,nt,nc,ne,0 if q else d-1,-be,-a,q)
  if sc>=be:return be
  a=max(a,sc)
 return a
def F(s):
 q=s.split();b=['.']*64;i=0
 for c in q[0]:
  if c=='/':continue
  elif c.isdigit():i+=int(c)
  else:b[i]=c;i+=1
 return b,(q[1]if len(q)>1 else'w'),(q[2]if len(q)>2 else''),(None if len(q)<4 or q[3]=='-'else(8-int(q[3][1]))*8+ord(q[3][0])-97)
def X(m):return chr(97+m[0]%8)+str(8-m[0]//8)+chr(97+m[1]%8)+str(8-m[1]//8)+(m[2].lower()if m[2]else'')
def B(b,tn,cas,ep,d=3):
 ms=L(b,tn,cas,ep);bv=-10**9;bm=ms[0]
 for m in ms:
  nb,nt,nc,ne=M(b,tn,cas,ep,m);sc=-S(nb,nt,nc,ne,d-1,-10**9,10**9)
  if sc>bv:bv,bm=sc,m
 return bm
Z='rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq -'
def main():
 b,tn,cas,ep=F(Z)
 for ln in sys.stdin:
  p=ln.split()
  if not p:continue
  c=p[0]
  if c=='uci':print('id name MicroSun\nuciok')
  elif c=='isready':print('readyok')
  elif c=='quit':break
  elif c=='position':
   if 'fen'in p:
    j=p.index('fen');k=p.index('moves')if 'moves'in p else len(p);b,tn,cas,ep=F(' '.join(p[j+1:k]))
   else:b,tn,cas,ep=F(Z)
   if 'moves'in p:
    for s in p[p.index('moves')+1:]:
     f=(8-int(s[1]))*8+ord(s[0])-97;t=(8-int(s[3]))*8+ord(s[2])-97;pr=s[4]if len(s)>4 else''
     m=next((x for x in L(b,tn,cas,ep)if x[0]==f and x[1]==t and x[2].lower()==pr),None)
     if m:b,tn,cas,ep=M(b,tn,cas,ep,m)
  elif c=='go':print('bestmove',X(B(b,tn,cas,ep,min(int(p[p.index('depth')+1])if 'depth'in p else 3,4))))
  sys.stdout.flush()
if __name__=='__main__':main()
