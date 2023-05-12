outfl=open("stations.lst","w")
for line in open("stainfo.dat","r").readlines():
    sta = line.split()[0]
    lat = line.split()[1]
    lon = float(line.split()[2])
    if lon < 0:
        lon += 360
    outfl.write("%s\t%s\t%f\n" %(sta,lat,lon) )

outfl.close()