#! /usr/bin/env python
import sys
import ROOT as rt

# read the starlight signal and background data and produced mixed data with specified signal to background ratio
if __name__ == "__main__":
    # if the user includes the flag -h or --help print the options
    if '-h' in sys.argv or '--help' in sys.argv:
        print ("Usage: %s [-input1 <starlight_signal.dat>] [-input2 <starlight_background.dat>] [-ratio <signal_to_background_ratio>] [-output <output_file_name> <dooutput 0 or 1>]" % sys.argv[0])
        sys.exit(1)
    #default values
    signal_file = "starlight_signal.dat"
    bkg_file = "starlight_bkg.dat"
    ratio = 2.
    outfile = "mixeddata.dat"
    dooutput =0 # decision to do output

    if "-input1" in sys.argv:
        p = sys.argv.index("-input1")
        signal_file = sys.argv[p+1]

    if "-input2" in sys.argv:
        p = sys.argv.index("-input2")
        bkg_file = sys.argv[p+1]
    if "-ratio" in sys.argv:
        p = sys.argv.index("-ratio")
        ratio = sys.argv[p+1]
    if "-output" in sys.argv:
        p = sys.argv.index("-output")
        outfile = sys.argv[p+1]
        dooutput = int(sys.argv[p+2])

    #particle_mass =
    print(signal_file)
    with open(signal_file) as ifile,open(bkg_file) as ifile2:
        #print("signal_file")
        Nsignal = 0
        Nbkg =  0
        Ngamma = 0
        signaltobkgfile = open("ratio.txt", 'w')
        mixfile = open(outfile, 'w')



        #sigmasslist=[]
        Mass_high = 0. # Higest mass of signal sample
        Mass_low =100.   # lowest mass of singal sample
        for line in ifile:
             print (line)
             p1 = rt.TLorentzVector() #particle one
             p2 = rt.TLorentzVector() #particle two
             #P = rt.TLorentzVector() # parent p1+p2
             lineVals = line.split()

             #lineVals = float(lineVals)
             p1.SetPxPyPzE(float(lineVals[0]),float(lineVals[1]),float(lineVals[2]),float(lineVals[3]))
             p2.SetPxPyPzE(float(lineVals[4]),float(lineVals[5]),float(lineVals[6]),float(lineVals[7]))
             if (dooutput==1):
                 mixfile.write(lineVals[0]+" "+lineVals[1]+" "+lineVals[2]+" "+lineVals[3]+" ")
                 mixfile.write(lineVals[4]+" "+lineVals[5]+" "+lineVals[6]+" "+lineVals[7]+" \n")
             P = p1+p2
             if Mass_high<P.M():
                 Mass_high = P.M()
             if Mass_low>P.M():
                 Mass_low = P.M()
             Nsignal = Nsignal+1
             #sigmasslist.append[P.M()]
        print ("high and low mass limit",Mass_high, Mass_low)
        for line in ifile2:
            lineVals = line.split()
            p1 = rt.TLorentzVector() #particle one
            p2 = rt.TLorentzVector() #particle two
            p1.SetPxPyPzE(float(lineVals[0]),float(lineVals[1]),float(lineVals[2]),float(lineVals[3]))
            p2.SetPxPyPzE(float(lineVals[4]),float(lineVals[5]),float(lineVals[6]),float(lineVals[7]))
            if (dooutput==1):
                mixfile.write(lineVals[0]+" "+lineVals[1]+" "+lineVals[2]+" "+lineVals[3]+" ")
                mixfile.write(lineVals[4]+" "+lineVals[5]+" "+lineVals[6]+" "+lineVals[7]+" \n")
            P = p1+p2
            Ngamma = Ngamma +1
            if(P.M()<=(Mass_high+0.06) and P.M()>= (Mass_low-0.06)):
                #print("passed mass is :", P.M())
                Nbkg = Nbkg+1
        calcratio =  2#float(Nsignal)/float(Nbkg)



        if (abs(calcratio-float(ratio))<1.):
            dooutput = 1
            print ("ratio of singal to background is : ", calcratio,"produce is: ", Nsignal, Nbkg, Ngamma)
            print("required ratio is:",ratio,"saving the files to ", outfile)
        else:
            dooutput = 0


        if (dooutput==0):

            print ("ratio of singal to background is : ", calcratio,"produce is: ", Nsignal, Nbkg, Ngamma)
            print("required ratio is:",ratio,"running process again")
        #print ("the output is ", dooutput)
        signaltobkgfile.write(str(dooutput)+" \n")
