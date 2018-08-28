set deviceNum=20
set sampleNum=5

for /L %%j in (1,1,%sampleNum%) do bin\Debug\DMaster ..\Sample\sample%deviceNum%\sample%deviceNum%_%%j.in  Result_%deviceNum%.out
pause