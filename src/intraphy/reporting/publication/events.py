"""Elementary model edits on retained or source-supported DNA."""
from .drawing import box, line, text, track, document, RULE, PALE, GOLD, TEAL, PURPLE, MUTED


def track_event(x, y, title, before, after, annotation):
    """A before/after glyph for one admissible geometry transition."""
    b = [box(x,y,530,242,"white",RULE), text(x+20,y+35,title,18,bold=True),
         text(x+18,y+105,"Before",12,color=MUTED), text(x+18,y+168,"After",12,color=MUTED)]
    gx=x+105; width=390
    b += track(gx,y+99,width,100,before,height=21)
    b += [line(gx+width/2,y+118,gx+width/2,y+143,GOLD,2,arrow=True)]
    b += track(gx,y+163,width,100,after,height=21)
    b += [text(x+20,y+218,annotation,12,color=MUTED)]
    return b


def draw():
    b = [text(42, 46, "Elementary model edits", 27, bold=True),
         text(42, 80, "Allowed transitions in the declared genomic exon-span edit graph", 14, color=MUTED)]
    b += track_event(42,105,"A   Split",[(8,92)],[(8,46),(54,92)],"Separator added; DNA retained")
    b += track_event(625,105,"B   Fusion",[(8,43),(57,92)],[(8,92)],"Intervening DNA retained")
    b += track_event(1208,105,"C   Boundary move",[(15,77)],[(25,77)],"One supported span edge moves")
    b += track_event(42,370,"D   Exonization",[(7,28),(72,93)],[(7,28),(42,58),(72,93)],"Candidate interval; DNA retained")
    b += track_event(625,370,"E   Exon inactivation",[(7,28),(42,58),(72,93)],[(7,28),(72,93)],"Annotation-conditional change")
    # The outlined interval is source-supported DNA entering the local tract.
    b += [box(1208,370,530,242,"white",RULE), text(1228,405,"F   DNA insertion",18,bold=True),
          text(1226,473,"Before",12,color=MUTED), text(1226,536,"After",12,color=MUTED)]
    gx=1305; width=390
    b += track(gx,464,width,100,[(8,42),(58,92)],height=21)
    b += [f'<rect x="{gx+width*.42}" y="461" width="{width*.16}" height="6" fill="white"/>',
          line(gx+width/2,483,gx+width/2,508,GOLD,2,arrow=True)]
    b += track(gx,528,width,100,[(8,42),(48,53),(58,92)],height=21)
    b += [f'<rect x="{gx+width*.42}" y="517" width="{width*.16}" height="28" fill="none" stroke="{GOLD}" stroke-width="2"/>',
          text(1228,598,"Material m: absent → present",12,color=MUTED)]
    # One before/after geometry: a single deleted DNA tract changes multiple
    # exon spans in that geometry; no coexisting transcript structures shown.
    b += [box(42,650,1696,300,PALE,RULE),
          text(64,690,"G   One DNA deletion changes multiple spans in one geometry",17,bold=True),
          text(72,735,"Before",12,color=MUTED), text(838,735,"After",12,color=MUTED)]
    left_x,left_w=135,535
    source_x=left_x+left_w*.36; source_w=left_w*.24
    b += [f'<rect x="{source_x}" y="754" width="{source_w}" height="36" fill="{GOLD}" fill-opacity=".22"/>']
    b += track(left_x,766,left_w,100,[(4,30),(36,52),(56,72),(80,96)],height=21)
    b += [text(left_x,810,"One geometry",12,color=TEAL,bold=True),
          line(702,766,790,766,GOLD,3,arrow=True)]
    right_x,right_w=885,500
    b += track(right_x,766,right_w,100,[(4,30),(60,72),(80,96)],height=21)
    # Mask the deleted segment on the after DNA backbone.
    gap_x=right_x+right_w*.36; gap_w=right_w*.24
    b += [f'<rect x="{gap_x}" y="763" width="{gap_w}" height="6" fill="{PALE}"/>',
          text(right_x,810,"Changed exon-span geometry",12,color=PURPLE,bold=True),
          text(1450,790,"Endpoint contrast",13,bold=True),
          text(1450,829,"does not specify",12,color=MUTED),
          text(1450,863,"path or mechanism",12,color=MUTED)]
    b += [text(64,918,"Edit motifs are model transitions conditional on the declared catalogue and annotation support.",12,color=MUTED)]
    return document(1790, 985, "Elementary model edits", b,
                    "Panels A–F illustrate allowed model transitions including split and fusion on retained DNA, exon-boundary movement, exonization, annotation-conditional exon inactivation, and source-supported DNA insertion. Panel G shows one before/after genomic exon geometry in which a single DNA deletion changes multiple exon spans; it does not show coexisting transcript configurations. The motifs are conditional on the declared catalogue and annotation support. Branch endpoint labels do not identify a transition path, event count, or molecular mechanism.")
