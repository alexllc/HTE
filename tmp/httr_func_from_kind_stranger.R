library(httr)  #Talks to the web server

#Call the web server to list all pathways and store in resp
resp <- GET("http://rest.kegg.jp/list/pathway")

#resp is a tab separated value (tsv). Kind of annoying, since everyone expects JSON strings
#from webservers. So we have to do a bit of an odd thing

#R doesn't have an easy function to read a tsv from a stored string
#So instead, save the file down
write(content(resp, "text"), "out.tsv")

#Then load the file in using core R. No headers.
tDat = read.delim("out.tsv", header=FALSE, sep="\t")

#Your data is now in the dataframe 'tDat' with V1 being the path and V2 being the description