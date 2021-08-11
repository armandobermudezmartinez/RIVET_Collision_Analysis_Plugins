BEGIN PLOT /CMS_2020_I1776758/d*
#Title=[Insert title for histogram d01-x01-y01 here]
# + any additional plot settings you might like, see make-plots documentation
END PLOT

BEGIN PLOT /CMS_2020_I1776758/d01-x01-y01
Title=CMS - Ratio charm/jet (combined) 
YLabel=$R(c/j)$
XLabel=$p_T (jet)$ [GeV]
YMax=0.16
YMin=0.0
RatioPlotYMax=2.0
RatioPlotYMin=0.0
LogX=0
LogY=0
# + any additional plot settings you might like, see make-plots documentation
END PLOT
BEGIN PLOT /CMS_2020_I1776758/d03-x01-y01
Title=CMS - Ratio bottom/jet (combined) 
YLabel=$R(b/j)$
XLabel=$p_T (jet)$ [GeV]
YMax=0.125
YMin=0.0
RatioPlotYMax=2.0
RatioPlotYMin=0.0
LogX=0
LogY=0
# + any additional plot settings you might like, see make-plots documentation
END PLOT
BEGIN PLOT /CMS_2020_I1776758/d05-x01-y01
Title=CMS - Ratio charm/bottom (combined) 
YLabel=$R(c/b)$
XLabel=$p_T (jet)$ [GeV]
YMax=2.6
YMin=0.0
RatioPlotYMax=2.0
RatioPlotYMin=0.0
LogX=0
LogY=0
# + any additional plot settings you might like, see make-plots documentation
END PLOT


BEGIN PLOT /CMS_2020_I1776758/d02-x01-y01
Title=CMS - Ratio charm/jet (combined)
YLabel=$R(c/j)$
XLabel=$p_T (Z)$ [GeV] 
YMax=0.16
YMin=0.0
RatioPlotYMax=2.0
RatioPlotYMin=0.0
LogX=0
LogY=0
# + any additional plot settings you might like, see make-plots documentation
END PLOT
BEGIN PLOT /CMS_2020_I1776758/d04-x01-y01
Title=CMS - Ratio bottom/jet (combined)
YLabel=$R(b/j)$
XLabel=$p_T (Z)$ [GeV] 
YMax=0.125
YMin=0.0
RatioPlotYMax=2.0
RatioPlotYMin=0.0
LogX=0
LogY=0
# + any additional plot settings you might like, see make-plots documentation
END PLOT
BEGIN PLOT /CMS_2020_I1776758/d06-x01-y01
Title=CMS - Ratio charm/bottom (combined)
YLabel=$R(c/b)$
XLabel=$p_T (Z)$ [GeV] 
YMax=2.6
YMin=0.0
RatioPlotYMax=2.0
RatioPlotYMin=0.0
LogX=0
LogY=0
# + any additional plot settings you might like, see make-plots documentation
END PLOT

# ... add more histograms as you need them ...
