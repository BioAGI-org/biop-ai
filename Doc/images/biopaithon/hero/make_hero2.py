"""BioPAIthon hero, take two: a thin double helix that becomes a circuit as it crosses the frame."""
import math
W,H=1600,900; CY=450; A=165; LAM=420; X0=60; X1=1540
BG="#0b1424"; GOLD="#f4c542"; BLUE="#5aa9ef"; CYAN="#8fdcff"; DIM=0.55
SEAM=900                            # hand-off at a crossing of the two strands                       # transition zone: organic -> silicon
def sine(x,phi): return CY + A*math.sin(2*math.pi*(x-X0)/LAM + phi)
def depth(x,phi): return math.cos(2*math.pi*(x-X0)/LAM + phi)   # >0 means in front
def smooth(phi, xa, xb, step=8):
    xs=list(range(int(xa),int(xb)+1,step))+[xb]
    return "M"+" L".join(f"{x:.1f},{sine(x,phi):.1f}" for x in xs)
def traced(phi, xa, xb, step=16):
    """Manhattan-plus-45-degree polyline that tracks the sine within half a step."""
    x=xa; y=round(sine(xa,phi)); pts=[(x,y)]
    sign=lambda v: 1 if v>0 else -1
    while x<xb:
        target=sine(x+step,phi); diff=target-y
        while abs(diff)>step:                 # steep part: vertical run first
            y+=step*sign(diff); pts.append((x,y)); diff=target-y
        dy=step*sign(diff) if abs(diff)>step/2 else 0
        x+=step; y+=dy; pts.append((x,y))
    turns=[]; last=(-999,-999)
    for i in range(1,len(pts)-1):
        if (pts[i][0]-pts[i-1][0],pts[i][1]-pts[i-1][1])!=(pts[i+1][0]-pts[i][0],pts[i+1][1]-pts[i][1]):
            if abs(pts[i][0]-last[0])+abs(pts[i][1]-last[1])>=44:   # keep vias sparse
                turns.append(pts[i]); last=pts[i]
    return "M"+" L".join(f"{p[0]:.1f},{p[1]:.1f}" for p in pts), turns, pts
o=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">','<defs>',
 f'<radialGradient id="bg" cx="0.5" cy="0.5" r="0.75"><stop offset="0" stop-color="#15284a"/><stop offset="1" stop-color="{BG}"/></radialGradient>',
 f'<clipPath id="cOrganic"><rect x="0" y="0" width="{SEAM}" height="{H}"/></clipPath>',
 f'<clipPath id="cSilicon"><rect x="{SEAM}" y="0" width="{W-SEAM}" height="{H}"/></clipPath>',
 f'<radialGradient id="glow"><stop offset="0" stop-color="{CYAN}" stop-opacity="0.9"/><stop offset="0.4" stop-color="{BLUE}" stop-opacity="0.35"/><stop offset="1" stop-color="{BLUE}" stop-opacity="0"/></radialGradient>',
 f'<pattern id="dots" width="40" height="40" patternUnits="userSpaceOnUse"><circle cx="20" cy="20" r="1.3" fill="#7fa6d8" fill-opacity="0.35"/></pattern>',
 '</defs>', f'<rect width="{W}" height="{H}" fill="url(#bg)"/>',
 f'<rect x="{SEAM}" y="0" width="{W-SEAM}" height="{H}" fill="url(#dots)" clip-path="url(#cSilicon)"/>']
strands=[(0.0,GOLD),(math.pi,BLUE)]
# ---- rungs (base pairs): organic = soft lines with end dots; silicon = pin-like with square pads
rung_x=list(range(X0+30,X1-20,42))
def rungs(kind):
    s=[]
    for x in rung_x:
        y1,y2=sine(x,0.0),sine(x,math.pi)
        if abs(y1-y2)<34: continue
        if kind=="organic":
            s.append(f'<line x1="{x}" y1="{y1:.1f}" x2="{x}" y2="{y2:.1f}" stroke="#c7d7ee" stroke-opacity="0.28" stroke-width="2"/>')
            for y,c in ((y1,GOLD),(y2,BLUE)): s.append(f'<circle cx="{x}" cy="{y:.1f}" r="3.2" fill="{c}" fill-opacity="0.9"/>')
        else:
            ym=(y1+y2)/2
            s.append(f'<line x1="{x}" y1="{y1:.1f}" x2="{x}" y2="{y2:.1f}" stroke="#9db8dd" stroke-opacity="0.35" stroke-width="2"/>')
            s.append(f'<rect x="{x-5}" y="{ym-9:.1f}" width="10" height="18" rx="1.5" fill="#1b2f52" stroke="#9db8dd" stroke-opacity="0.7" stroke-width="1.5"/>')
            for y,c in ((y1,GOLD),(y2,BLUE)): s.append(f'<rect x="{x-4}" y="{y-4:.1f}" width="8" height="8" fill="{c}"/>')
    return "".join(s)
o.append(f'<g clip-path="url(#cOrganic)">{rungs("organic")}</g>')
o.append(f'<g clip-path="url(#cSilicon)">{rungs("silicon")}</g>')
# ---- strands: back strand first, then front segments with a background halo (the classic helix gap)
def strand_layers(kind):
    out=[]
    for pass_ in ("back","front"):
        for phi,col in strands:
            if kind=="organic":
                d,turns=smooth(phi,X0,X1),[]
            else:
                d,turns,_=traced(phi,X0,X1)
            # split by depth into intervals
            ivs=[]; x=X0; cur=None
            while x<=X1:
                f=(depth(x,phi)>0)==(pass_=="front")
                if f and cur is None: cur=x
                if not f and cur is not None: ivs.append((cur,x)); cur=None
                x+=2
            if cur is not None: ivs.append((cur,X1))
            op=1.0 if pass_=="front" else DIM
            for xa,xb in ivs:
                cid=f"k{kind[0]}{pass_[0]}{int(phi*100)}{int(xa)}"
                out.append(f'<clipPath id="{cid}"><rect x="{xa}" y="0" width="{xb-xa}" height="{H}"/></clipPath>')
                g=f'<g clip-path="url(#{cid})" opacity="{op}">'
                if pass_=="front": g+=f'<path d="{d}" fill="none" stroke="{BG}" stroke-width="16" stroke-linejoin="round"/>'
                if kind=="organic":
                    g+=f'<path d="{d}" fill="none" stroke="{col}" stroke-width="5" stroke-linecap="round"/>'
                else:
                    g+=f'<path d="{d}" fill="none" stroke="{col}" stroke-width="4" stroke-linejoin="miter"/>'
                    g+="".join(f'<circle cx="{tx}" cy="{ty}" r="5" fill="{col}"/><circle cx="{tx}" cy="{ty}" r="2" fill="{BG}"/>' for tx,ty in turns if xa<=tx<=xb)
                out.append(g+'</g>')
    return "".join(out)
o.append(f'<g clip-path="url(#cOrganic)">{strand_layers("organic")}</g>')
o.append(f'<g clip-path="url(#cSilicon)">{strand_layers("silicon")}</g>')
# ---- a few glowing vias on the silicon side, and soft glow at the organic end
import random; random.seed(3)
_,_,ptsB=traced(math.pi,X0,X1); _,_,ptsG=traced(0.0,X0,X1)
cand=[p for p in ptsB+ptsG if p[0]>SEAM+60]
for x,y in random.sample(cand,7):
    o.append(f'<circle cx="{x}" cy="{y}" r="22" fill="url(#glow)"/><circle cx="{x}" cy="{y}" r="3" fill="{CYAN}"/>')
o.append('</svg>')
open('hero2.svg','w').write("\n".join(o)); print("written")
