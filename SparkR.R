installed.packages()
rattle()
sc

head(mtcars)
"head(mtcars)"
"SparkRSQL.init(sc)"
sqlContext<-sparkRSQL.init(sc)
head("mtc")

sdf<- createDataFrame(sqlContext, mtcars)
printSchema(sdf)

SparkR::head(sdf) #print content using SparkR
SparkR::head(select(sdf,sdf$mpg)) #Select a Column
SparkR::head(SparkR:: filter(sdf,sdf$mpg<18)) #filter based on a condition

sdf$wtTon<- sdf$wt * 0.45 # easy operation on a Spark data
SparkR::head(sdf)

SparkR::head(summarize(groupBy(sdf,sdf$cyl), wtavg= avg(sdf$wtTon)))

#Now Database by SparkR

registerTempTable(sdf,"cars")
highgearcars<-sql(sqlContext, "select gear from cars where cyl >=4 and cyl<=9")
SparkR::head(highgearcars)




