getHGMDData <- function(geneName,
                        user='hgmd',
                        password='hgmd',
                        host='192.168.99.100',
                        port='3306',
                        dbname='hgmd') {
  library("DBI")
  library("RMySQL")
  # load also HGMD
  con=dbConnect(MySQL(),username=user,password=password,host=host,unix.socket=port,dbname=dbname)
  # get all mutations
  res <- dbSendQuery(con,sprintf("select * from mutnomen inner join allmut as h1 on h1.acc_num=mutnomen.acc_num where gene='%s'",geneName))
  rawdat<-fetch(res,n=-1)
  cleanHGVS <- paste0("c.",rawdat$hgvs[!is.na(rawdat$hgvs)])
  dbDisconnect(con)
  (list(hgvs=cleanHGVS,rawData=rawdat))
}
getHGMDData("TP53")
