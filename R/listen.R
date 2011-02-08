listen <-
function(init = TRUE){
	if(init){
		cat("Choose 'BREAK' in the File menu of the graphics window or press CTRL+SHIFT+B if you wish to return to the console.")
	}
while (!is.null(ievent.wait())) {
if (iset.sel.changed()) {
resort(listen = T)
return(invisible(TRUE))
}
}
}

