listen <-
function(){
while (!is.null(ievent.wait())) {
if (iset.sel.changed()) {
resort(listen = T)
return(0)
}
}
}

