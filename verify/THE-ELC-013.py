import sys; sys.path.insert(0, '/home/claude/chem/work/g05'); from g05lib import check_item
E={'Zn':-0.76,'Fe':-0.44,'Cu':0.34}
anode=lambda a,b: a if E[a]<E[b] else b
assert anode('Fe','Cu')=='Fe' and anode('Fe','Zn')=='Zn'
import json; d=json.load(open('/home/claude/kurowasan3/alchemia/content/chem/ions.json'))
assert d['results']['Fe^2+']['Fe(CN)6-3']['color']=='濃青'
assert not check_item('THE-ELC-013')
