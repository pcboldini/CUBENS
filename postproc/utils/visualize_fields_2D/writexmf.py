def writexmf(filename,precision,x,y,z,timestamp,dt,dataNames):

    import numpy as np
    
    imax = np.size(x)
    jmax = np.size(y)
    kmax = np.size(z)
    tmax = np.size(timestamp)

    prec = 8

    if precision == 'single':
        prec = 4 

    XMLItemDesc  = 'Format=\"XML\" DataType=\"Float\" Precision=\"{}\"'.format(prec)
    dataItemDesc = 'Format=\"Binary\" DataType=\"Float\" Precision=\"{}\" Endian=\"Native\"'.format(prec)

    f = open(filename, 'wt')

# header
    f.write("<?xml version=\"1.0\" ?>\n")
    f.write("<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n")
    f.write("<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n")

# domain
    f.write("<Domain>\n")

    f.write("    <Topology name=\"TOPO\" TopologyType=\"3DRectMesh\" Dimensions=\"{0:5d}{1:5d}{2:5d}\"/>\n".format(kmax,jmax,imax))

# x, y, z coordinates 
    f.write("    <Geometry name=\"GEO\" GeometryType=\"VXVYVZ\">\n")
    f.write("        <DataItem {0} Dimensions=\"{1:5d}\">\n".format(XMLItemDesc,imax))
    for v in x: f.write('%15.6lE'% v)
    f.write("\n        </DataItem>\n")
    f.write("        <DataItem {0} Dimensions=\"{1:5d}\">\n".format(XMLItemDesc,jmax))
    for v in y: f.write('%15.6lE'% v)
    f.write("\n        </DataItem>\n")
    f.write("        <DataItem {0} Dimensions=\"{1:5d}\">\n".format(XMLItemDesc,kmax))
    for v in z: f.write('%15.6lE'% v)
    f.write("\n        </DataItem>\n")
    f.write("    </Geometry>\n")

# time coordinates
    f.write("    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n")
    f.write("        <Time TimeType=\"List\">\n")
    f.write("            <DataItem {0} Dimensions=\"{1:5d}\">\n".format(XMLItemDesc,tmax))
    for v in timestamp: f.write('%15.6lE'% (v*dt))
    f.write("\n            </DataItem>\n")
    f.write("        </Time>\n")
    
# files
    for v in timestamp: 
        f.write("        <Grid Name=\"T{0:07d}\" GridType=\"Uniform\">\n".format(v))
        f.write("            <Topology Reference=\"/Xdmf/Domain/Topology[1]\"/>\n")
        f.write("            <Geometry Reference=\"/Xdmf/Domain/Geometry[1]\"/>\n")
        for n in dataNames:
            f.write("            <Attribute Name=\"{0}\" Center=\"Node\">\n".format(n))
            f.write("                <DataItem {0} Dimensions=\"{1:5d}{2:5d}{3:5d}\">\n".format(dataItemDesc,kmax,jmax,imax))
            f.write("                    {0}.{1:07d}.bin\n".format(n,v))
            f.write("                </DataItem>\n")
            f.write("            </Attribute>\n")
        f.write("        </Grid>\n")
    
    
    f.write("    </Grid>\n")
    f.write("</Domain>\n")
    f.write("</Xdmf>\n")

    f.close()


