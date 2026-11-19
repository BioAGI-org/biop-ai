"""Generate the BioPAIthon blog hero: two helix snakes, one scaled, one silicon."""
import math, random, re
W,H=1600,900
cy=450; A=170; lam=520; x0=100; x1=x0+int(1.5*lam)
Y_PHI=math.pi/2; B_PHI=3*math.pi/2          # start at the extremes, so the tails are apart
def pt(x,phi): return (x, cy + A*math.sin(2*math.pi*(x-x0)/lam + phi))
def path(phi, xa, xb, dy=0, step=6):
    xs=[xa+i*step for i in range(int((xb-xa)/step)+1)]
    if xs[-1]<xb: xs.append(xb)
    return "M"+" L".join(f"{x:.1f},{pt(x,phi)[1]+dy:.1f}" for x in xs)
def front_intervals(phi):
    ivs=[]; x=x0; cur=None
    while x<=x1:
        f=math.cos(2*math.pi*(x-x0)/lam+phi)>0
        if f and cur is None: cur=x
        if not f and cur is not None: ivs.append((cur,x)); cur=None
        x+=2
    if cur is not None: ivs.append((cur,x1))
    return ivs
def ribbon_yellow(xa,xb):
    return (f'<path d="{path(Y_PHI,xa,xb)}" stroke="#4a3608" stroke-width="86" fill="none" stroke-linecap="round"/>'
            f'<path d="{path(Y_PHI,xa,xb)}" stroke="#f2c230" stroke-width="74" fill="none" stroke-linecap="round"/>'
            f'<path d="{path(Y_PHI,xa,xb)}" stroke="url(#scales)" stroke-width="74" fill="none" stroke-linecap="round"/>'
            f'<path d="{path(Y_PHI,xa,xb,-22)}" stroke="#ffe27a" stroke-width="10" fill="none" stroke-linecap="round" stroke-opacity="0.55"/>'
            f'<path d="{path(Y_PHI,xa,xb,26)}" stroke="#b07f14" stroke-width="12" fill="none" stroke-linecap="round" stroke-opacity="0.45"/>')
def ribbon_blue(xa,xb):
    s=(f'<path d="{path(B_PHI,xa,xb)}" stroke="#050b16" stroke-width="86" fill="none" stroke-linecap="round"/>'
       f'<path d="{path(B_PHI,xa,xb)}" stroke="#15294a" stroke-width="74" fill="none" stroke-linecap="round"/>')
    for dy,wd,op in ((-20,2.5,0.95),(-7,2,0.7),(7,2,0.7),(20,2.5,0.95)):
        s+=f'<path d="{path(B_PHI,xa,xb,dy)}" stroke="#d98e3a" stroke-width="{wd}" fill="none" stroke-opacity="{op}" stroke-linecap="round"/>'
    x=xa+30
    while x<xb-20:
        y=pt(x,B_PHI)[1]
        for dy in (-20,20): s+=f'<circle cx="{x:.1f}" cy="{y+dy:.1f}" r="5" fill="#f2c230" stroke="#4a3608" stroke-width="1.5"/>'
        gy=pt(x+40,B_PHI)[1]; s+=f'<circle cx="{x+40:.1f}" cy="{gy:.1f}" r="16" fill="url(#glow)"/><circle cx="{x+40:.1f}" cy="{gy:.1f}" r="3.5" fill="#cfefff"/>'
        x+=80
    return s
def head(phi, fill, outline, eye):
    x,y=pt(x1,phi); dy=2*math.pi*A/lam*math.cos(2*math.pi*(x1-x0)/lam+phi); ang=math.degrees(math.atan2(dy,1))
    return (f'<g transform="translate({x:.1f},{y:.1f}) rotate({ang:.1f})">'
            f'<path d="M-34,-40 C10,-52 70,-40 96,-12 C102,-5 102,5 96,12 C70,40 10,52 -34,40 Z" fill="{fill}" stroke="{outline}" stroke-width="6" stroke-linejoin="round"/>'
            f'<ellipse cx="40" cy="-17" rx="6" ry="8" fill="{eye}"/><ellipse cx="40" cy="17" rx="6" ry="8" fill="{eye}"/>'
            f'<path d="M100,0 l30,-9 M100,0 l30,9" stroke="#e2574c" stroke-width="4" fill="none" stroke-linecap="round"/></g>')
o=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">']
clips=[]; groups=[]
groups.append(ribbon_yellow(x0,x1)); groups.append(ribbon_blue(x0,x1))
for phi,fn in ((Y_PHI,ribbon_yellow),(B_PHI,ribbon_blue)):
    for xa,xb in front_intervals(phi):
        i=len(clips); clips.append(f'<clipPath id="c{i}"><rect x="{xa}" y="0" width="{xb-xa}" height="{H}"/></clipPath>')
        groups.append(f'<g clip-path="url(#c{i})">{fn(max(x0,xa-60),min(x1,xb+60))}</g>')
o.append('<defs><radialGradient id="glow"><stop offset="0" stop-color="#9fd8ff" stop-opacity="0.95"/><stop offset="0.5" stop-color="#4fa8e8" stop-opacity="0.45"/><stop offset="1" stop-color="#4fa8e8" stop-opacity="0"/></radialGradient>'
 '<radialGradient id="vig" cx="0.35" cy="0.5" r="0.8"><stop offset="0" stop-color="#132742"/><stop offset="1" stop-color="#070e1b"/></radialGradient>'
 '<pattern id="scales" width="26" height="20" patternUnits="userSpaceOnUse"><path d="M0,20 a13,13 0 0 1 26,0 M-13,10 a13,13 0 0 1 26,0 M13,10 a13,13 0 0 1 26,0" fill="none" stroke="#c99a1e" stroke-width="2.2" stroke-opacity="0.75"/></pattern>'
 '<pattern id="grid" width="48" height="48" patternUnits="userSpaceOnUse"><path d="M48,0 L0,0 0,48" fill="none" stroke="#1a2f4f" stroke-width="1"/></pattern>'+"".join(clips)+'</defs>')
o.append(f'<rect width="{W}" height="{H}" fill="url(#vig)"/><rect width="{W}" height="{H}" fill="url(#grid)" opacity="0.55"/>')
random.seed(7); pads=[]; o.append('<g stroke="#223f6b" stroke-width="2" fill="none" opacity="0.55">')
for i in range(22):
    x=random.randint(40,W-40); y=random.randint(40,H-40); pts=[(x,y)]
    for k in range(random.randint(2,4)):
        d=random.choice([(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,1),(1,-1),(-1,-1)]); L=random.randint(60,220); x,y=x+d[0]*L,y+d[1]*L; pts.append((x,y))
    o.append('<polyline points="'+" ".join(f"{p[0]},{p[1]}" for p in pts)+'"/>'); pads+=[pts[0],pts[-1]]
o.append('</g><g fill="#0b1526" stroke="#2b4f82" stroke-width="2" opacity="0.7">'+"".join(f'<circle cx="{x}" cy="{y}" r="5"/>' for x,y in pads)+'</g>')
o+= [f'<g>{g}</g>' if not g.startswith('<g') else g for g in groups]
o.append(head(Y_PHI,"#f2c230","#4a3608","#1b1b1b")); o.append(head(B_PHI,"#15294a","#050b16","#8fd3ff")); o.append('</svg>')
open('hero_vector.svg','w').write("\n".join(o)); print("svg written")
