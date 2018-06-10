bhts2HTML = function(dat, dir, fname, title=NULL, bgcolor="#BBBBEE"){
  dat = data.frame(compound=rownames(dat), dat)
  rownames(dat) = NULL
  target = HTMLInitFile(dir, filename=fname, Title=title, BackGroundColor=bgcolor)
  HTML(xtable(dat), file=target)
  HTMLEndFile()
}
        
