##!/usr/bin/env python
import numpy as np
from scipy.constants import c, pi, e, m_e, epsilon_0
from math import ceil
from datetime import date, datetime
import sys
E0 = 0.5110  # [MeV] electron rest energy
E0p = 938.08 # [MeV] proton rest energy

def read_default_parameters(input_file):
    parameters = {}
    with open(input_file, "r") as file:
        for line in file:
            # Ignore lines starting with #
            if line.startswith("#") or line.startswith("&"):
                continue

            # Split line by '#' to ignore comments
            line_parts = line.split("#", 1)
            line = line_parts[0].strip()  # Take the first part before the '#'


            # Split the line into key and value
            key_value = line.strip().split("=")

            # Ensure the line has the correct format (contains key and value)
            if len(key_value) != 2:
                continue

            # Extract key and value
            key, value = key_value

            # Add key-value pair to parameters dictionary
            parameters[key.strip()] = value.strip()

    # Convert values to the right variable
    for key, value in parameters.items():

        try:
            # First, try converting to an integer
            result = int(value)
            parameters[key] = result #, "integer"
        except ValueError:
            try:
                # If it fails, try converting to a float
                result = float(value)
                parameters[key] = result #, "float"
            except ValueError:
                # If both fail, it remains a string
                parameters[key] = value #, "string"

            #############
    return parameters


def generat_inputdeck(parameters, input_deck):
    with open(input_deck, "w") as file:

        comment_plasma = "# Plasma: "
        comment_driver = "# Driver: "
        comment_witness= "# Witness: "
        comment_radiation = "# Betatron: "
        comment_simulation = "# Simulation: "
        today = date.today()
        d1 = today.strftime("%d/%m/%Y")
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
        About_inputs = parameters.get("About_inputs")
        file.write("# ************************************************\n")
        file.write("# Automatic QV3D input deck\n")
        file.write("# Date: " + dt_string + "\n")
        file.write("# About: " + About_inputs + "\n")
        file.write("# ************************************************\n\n")
        file.write("&Domain\n")
        # #####################
        file.write("  # Moving window length in normalized units (simulation domain)\n")
        Xlength = parameters.get("Xlength")
        file.write("  Xlength = " + str(Xlength) + "\n")
        Ylength = parameters.get("Ylength")
        file.write("  Ylength = " + str(Ylength) + "\n")
        Zlength = parameters.get("Zlength")
        file.write("  Zlength = " + str(Zlength) + "\n")
        comment_simulation += "(x,y,z)=("+ str("{:.0f}".format(Xlength))+ ','+ str("{:.0f}".format(Ylength))+ ','+ str("{:.0f}".format(Zlength))+'); '

        file.write("\n  # Length resolution - cell size in normalized units\n")
        Hx = parameters.get("Hx")
        file.write("  Hx = " + str(Hx) + "\n")
        NxSplit = parameters.get("NxSplit")
        file.write("  NxSplit = " + str(int(NxSplit)) + "\n")
        Hy = parameters.get("Hy")
        file.write("  Hy = " + str(Hy) + "\n")
        Hz = parameters.get("Hz")
        file.write("  Hz = " + str(Hz) + "\n")
        comment_simulation += "(dx,dy,dz)=("+ str("{:.3f}".format(Hx))+ ','+ str("{:.3f}".format(Hy))+ ','+ str("{:.3f}".format(Hz))+'); '

        nIter = parameters.get("nIter")
        file.write("\n  nIter = " + str(int(nIter)) + "\n")
        Accuracy = parameters.get("Accuracy")
        file.write("  Accuracy = " + str(Accuracy) + "\n")

        file.write("\n  # Time parameters\n")
        Ts = parameters.get("Ts")
        file.write("  Ts = " + str(Ts) + "			# This is the time step.\n")
        comment_simulation += "dt=" + str(Ts) + '; '
        TsIni = parameters.get("TsIni")
        file.write("  TsIni = " + str(TsIni) + "\n")
        TimeIni = parameters.get("TimeIni")
        file.write("  TimeIni = " + str(TimeIni) + "\n\n")


        file.write("\n  # Particles\n")
        n0 = parameters.get("plasma_density")
        n0*=1e6  	# convert density to /m3
        comment_plasma += "n0=" + str("{:.2e}".format(n0/1e6))+"/cm3; "
        wp = np.sqrt(n0*e**2/epsilon_0/m_e)
        kp = wp/c
        lamdap = 100*2*pi*c/wp 			# plasma wavelength in cm
        file.write("  Wavelength = " + str("{:.4e}".format(lamdap)) + "	# plasma wavelength [SHOULD be in cm]\n")

        Nspecies = parameters.get("Nspecies")
        file.write("  Nspecies = " + str(int(Nspecies)) + "			# includes plasma, driver, witness\n")
        Npulses = parameters.get("Npulses")
        file.write("  Npulses = " + str(int(Npulses)) + "\n")

        file.write("\n  # Output diag\n")
        NMovieFramesH5 = parameters.get("NMovieFramesH5")
        file.write("  NMovieFramesH5 = " + str(int(NMovieFramesH5)) + "		# Output in 3d: 0 variables (see &MovieHDF5)\n")
        NMovie2dFramesH5 = parameters.get("NMovie2dFramesH5")
        file.write("  NMovie2dFramesH5 = " + str(int(NMovie2dFramesH5)) + "		# Output in 2d: 6 variables (see &Movie2dHDF5)\n")

        file.write("\n  # PIC algorithm\n")
        MaxwellSolver = parameters.get("MaxwellSolver")
        file.write("  MaxwellSolver = " + str(int(MaxwellSolver)) + "		# 0 -> RFFT; 1 -> CFFT\n")
        ParticlePusher = parameters.get("ParticlePusher")
        file.write("  ParticlePusher = " + str(int(ParticlePusher)) + "		# 0 -> NGP, 1-> trilinear, 2-> quadratic, 15 -> quadratic interpolation in the transverse direction \n")
        FollowParticles = parameters.get("FollowParticles")
        file.write("  FollowParticles = " + str(int(FollowParticles)) + "		# in order to save plasma particles; advance topic - better not use it (0 to deactivate, 1 to activate) \n")
        Guess = parameters.get("Guess")
        file.write("  Guess = " + str(int(Guess)) + "			# 0 - simple, 1 - linear\n")
        Hybrid = parameters.get("Hybrid")
        file.write("  Hybrid = " + str(int(Hybrid)) + "			# 1-linear; 2-density change; 3-momenta nonlinearity and B-response\n")
        DiffractLasers = parameters.get("DiffractLasers")
        file.write("  DiffractLasers = " + str(int(DiffractLasers)) + "\n")
        RefractLasers = parameters.get("RefractLasers")
        file.write("  RefractLasers = " + str(int(RefractLasers)) + "\n")
        NSplitTimeStepRear = parameters.get("NSplitTimeStepRear")
        file.write("  NSplitTimeStepRear = " + str(int(NSplitTimeStepRear)) + "\n")
        RandomSeed = parameters.get("RandomSeed")
        file.write("  RandomSeed = " + str(int(RandomSeed)) + "		# 0: grid loading - 1: particles initialised with random position inside cell\n")
        RearPosition = parameters.get("RearPosition")
        file.write("  RearPosition = " + str(int(RearPosition)) + "		# 3.5 advance topic- to make use of ContinueBack otherwise should be set to a negative value \n")
        SuppressLogs = parameters.get("SuppressLogs")
        file.write("  SuppressLogs = " + str(int(SuppressLogs)) + "\n")
        file.write("/\n# ************************************************\n")

        file.write("&MPP_partition\n")
        #######################
        file.write("  # Maximum should be along X\n")
        file.write("  # The product is number of cores you choose\n")
        file.write("  # example for toy model 100, 2, 2\n")

        Xpartition = parameters.get("Xpartition")
        Ypartition = parameters.get("Ypartition")
        Zpartition = parameters.get("Zpartition")
        file.write("  Xpartition = " + str(int(Xpartition)) + "\n")
        file.write("  Ypartition = " + str(int(Ypartition)) + "\n")
        file.write("  Zpartition = " + str(int(Zpartition)) + "\n")
        comment_simulation += "CPUs="+str(int(Xpartition*Ypartition*Zpartition))+"=("+str(int(Xpartition))+','+str(int(Ypartition))+','+str(int(Zpartition))+'); '
        file.write("/\n# ************************************************\n")

        file.write("&Controls\n")
        #######################
        ContinueBack = parameters.get("ContinueBack")
        file.write("  ContinueBack = " + str(int(ContinueBack)) + "		# advanced topic\n")
        Reload = parameters.get("Reload")
        file.write("  Reload = " + str(int(Reload)) + "			# to resume form input data\n")
        Nwrite = parameters.get("Nwrite")
        file.write("  Nwrite = " + str(int(Nwrite)) + "			# to resume form input data\n")
        FirstSaveTime = parameters.get("FirstSaveTime")
        file.write("  FirstSaveTime = " + str(int(FirstSaveTime)) + "		# to alter the default: first movie frame is saved after the first time step; first particle saves is after SavePeriod \n")
        PostProcessing = parameters.get("PostProcessing")
        file.write("  PostProcessing = " + str(int(PostProcessing)) + "		# code internal postprocessing plots vlp3\n")
        Ndiagnose = parameters.get("Ndiagnose")
        file.write("  Ndiagnose = " + str(int(Ndiagnose)) + "\n")
        CPUstop = parameters.get("CPUstop")
        file.write("  CPUstop = " + str(int(CPUstop)) + "\n")

        plasma_length = parameters.get("plasma_length")
        comment_plasma+= 'Length='+ str(plasma_length) + 'm; '
        output_every_what_distance = parameters.get("output_every_what_distance")
        comment_plasma+='diag every ' + str(output_every_what_distance) + 'm; '
        PhaseStop = plasma_length*wp/c
        number_of_outputs_along_plasma = int( round(plasma_length/output_every_what_distance) )
        SavePeriod = int( round(output_every_what_distance*wp/c) )
        PhaseStop_rounded = SavePeriod * number_of_outputs_along_plasma
        file.write("\n  # Run time and the outputs periods\n")
        file.write("  PhaseStop = " + str(PhaseStop_rounded) + "	 	# (" + str("{:.2f}".format(PhaseStop)) + "=PhaseStop=plasmaLength*wp/c). Rounded to have an int SavePeriod. (Rounded plasma length="+ str("{:.2f}".format(PhaseStop_rounded*c/wp)) +   " vs Real ="+ str("{:.2f}".format(plasma_length)) +" m.)\n")


        file.write("  SavePeriod = " + str(SavePeriod) + "		# Save particles in h5files/vs???_3d_particles.h5 files- every " + str("{:.1f}".format(output_every_what_distance*100)) + " cm inside plasma\n")
        Movie2dPeriodH5 = SavePeriod
        file.write("  Movie2dPeriodH5 = " + str(Movie2dPeriodH5) + "	# saves densities, fields etc in h5files/v2d_mframe_?????.h5 every " + str("{:.1f}".format(output_every_what_distance*100)) + " cm inside plasma\n")
        SaveFieldsFlag = parameters.get("SaveFieldsFlag")
        file.write("  SaveFieldsFlag = " + str(int(SaveFieldsFlag)) + "\n")
        FullSaveStride = parameters.get("FullSaveStride")
        file.write("  FullSaveStride = " + str(int(FullSaveStride)) + "\n")
        file.write("/\n# ************************************************\n")

        file.write("&Pulse0\n")
        ######################
        a0 = parameters.get("a0")
        file.write("  a0 = " + str(a0) + "\n")
        Ypol = parameters.get("Ypol")
        file.write("  Ypol = " + str(Ypol) + "\n")
        Zpol = parameters.get("Zpol")
        file.write("  Zpol = " + str(Zpol) + "\n")
        Length = parameters.get("Length")
        file.write("  Length = " + str(int(Length)) + "\n")
        Ywidth = parameters.get("Ywidth")
        file.write("  Ywidth = " + str(int(Ywidth)) + "\n")
        Zwidth = parameters.get("Zwidth")
        file.write("  Zwidth = " + str(int(Zwidth)) + "\n")
        RiseTime = parameters.get("RiseTime")
        file.write("  RiseTime = " + str(int(RiseTime)) + "\n")
        DropTime = parameters.get("DropTime")
        file.write("  DropTime = " + str(int(DropTime)) + "\n")
        Xcenter = parameters.get("Xcenter")
        file.write("  Xcenter = " + str(int(Xcenter)) + "\n")
        Xperiod = parameters.get("XperiodPulse")
        file.write("  Xperiod = " + str(int(Xperiod)) + "\n")
        Ycenter = parameters.get("Ycenter")
        file.write("  Ycenter = " + str(int(Ycenter)) + "\n")
        Zcenter = parameters.get("Zcenter")
        file.write("  Zcenter = " + str(int(Zcenter)) + "\n")
        Yphase = parameters.get("Yphase")
        file.write("  Yphase = " + str(Yphase) + "\n")
        Zphase = parameters.get("Zphase")
        file.write("  Zphase = " + str(Zphase) + "\n")
        FromBoundary = parameters.get("FromBoundary")
        file.write("  FromBoundary = " + str(FromBoundary) + "\n")
        Kx = parameters.get("Kx")
        file.write("  Kx = " + str(int(Kx)) + "\n")
        Ky = parameters.get("Ky")
        file.write("  Ky = " + str(int(Ky)) + "\n")
        Kz = parameters.get("Kz")
        file.write("  Kz = " + str(int(Kz)) + "\n")
        file.write("/\n# ************************************************\n")

        file.write("&Electrons       	     	# This is always the background plasma.\n")
        ################################################################################
        Distribution = parameters.get("Distribution")
        file.write("  Distribution = " + str(int(Distribution)) + "		# (=9) uniform plasma [see plasma.cpp file in the source files]\n")
        if Distribution == 9:
            comment_plasma+= "uniform distribution; "
        else:
            raise ValueError("check plasma distribution. change plasma comment.")
        file.write("  Density = " + str(n0/n0) + "			# Normalized to the plasma density. n0=" + str("{:.2e}".format(n0/1e6)) + "/cm3\n")
        Begin = parameters.get("Begin")
        file.write("  Begin = " + str(Begin) + "			# Beginning of plasma on x- left boundary \n")
        PlateauBegin = parameters.get("PlateauBegin")
        file.write("  PlateauBegin = " + str(int(PlateauBegin)) + "		# beginning of polateau on x\n")
        PlateauEnd = parameters.get("PlateauEnd")
        file.write("  PlateauEnd = " + str("{:.1e}".format(PlateauEnd)) + "		# End of plateau on x \n")
        End = PhaseStop_rounded + 1000
        file.write("  End = " + str(int(End)) + "			# (PhaseStop="+str("{:.0f}".format(PhaseStop_rounded))+") end of plasma - right boundary at " + str("{:.1f}".format(plasma_length)) +" m.\n")
        Xperiod = parameters.get("Xperiod")
        file.write("  Xperiod = " + str(int(Xperiod)) + "\n")
        Delta = parameters.get("Delta")
        file.write("  Delta = " + str(int(Delta)) + "\n")
        RadiusY = Ylength
        RadiusZ = Zlength
        file.write("  RadiusY = " + str(RadiusY) + "	          	# = Ylength=" + str((Ylength)) +"\n")
        file.write("  RadiusZ = " + str(RadiusZ) + "	           	# = Zlength=" + str((Zlength)) +"\n")
        Px0 = parameters.get("Px0")
        Py0 = parameters.get("Py0")
        Pz0 = parameters.get("Pz0")
        file.write("  Px0 = " + str(Px0) + "			#initial momentum; cold plasma\n")
        file.write("  Py0 = " + str(Py0) + "			#initial momentum; cold plasma\n")
        file.write("  Pz0 = " + str(Pz0) + "			#initial momentum; cold plasma\n")
        P_perCell = parameters.get("P_perCell")
        comment_plasma += 'ppc='+str(int(P_perCell))
        file.write("  P_perCell = " + str(int(P_perCell)) + "			#particelpercell at the peak density (choose a square number)\n")

        file.write("/\n# ************************************************\n")

        file.write("&Specie1 			# Driver beam\n")
        #################################################
        Distribution = parameters.get("Distribution1")
        if Distribution == 1:
            comment_driver+= " Gaussian profile; "
        else:
            raise ValueError("check driver distribution. change driver comment.")
        file.write("  Distribution = " + str(int(Distribution)) + "		# Gaussian ellipsoid\n")
        Q = parameters.get("Charge1") * 1e-12 #C
        comment_driver+= 'q='+str("{:.2f}".format(Q/1e-12)) + 'pC; '
        sigx = parameters.get("sigx1")
        comment_driver+= 'sigx='+str("{:.2f}".format(sigx/1e-6)) + 'um; '
        sigr = parameters.get("sigr1")
        comment_driver+= 'sigr='+str("{:.2f}".format(sigr/1e-6)) + 'um; '
        ne = Q/(2*pi)**(3./2.)/sigx/sigr**2/e
        comment_driver+= 'n/n0='+str("{:.2f}".format(ne/n0)) + '; '
        file.write("  Density = " + str("{:.2e}".format(ne/n0)) + "		# normalized to the plasma density; ne=" + str("{:.2e}".format(ne/1e6)) +"/cm3 and n0=" + str("{:.2e}".format(n0/1e6)) +"/cm3\n")

        Begin = parameters.get("Begin1")
        file.write("  Begin = " + str(Begin) + "\n")
        PlateauBegin = parameters.get("PlateauBegin1")
        file.write("  PlateauBegin = " + str(int(PlateauBegin)) + "\n")
        PlateauEnd = parameters.get("PlateauEnd1")
        file.write("  PlateauEnd = " + str("{:.1e}".format(PlateauEnd)) + "\n")
        End = parameters.get("End1")
        file.write("  End = " + str("{:.1e}".format(End)) + "\n")

        sigy = sigz = sigr
        #print(sigx, kp)
        file.write("  RadiusX = " + str("{:.4f}".format( np.sqrt(2.)*sigx*kp )) + "		# =sqrt(2)sigx*kp - sigx = " + str("{:.2f}".format(sigx/1e-6)) + "um\n")
        file.write("  RadiusY = " + str("{:.4f}".format( np.sqrt(2.)*sigy*kp )) + "		# =sqrt(2)sigy*kp - sigy = sigr = " + str("{:.2f}".format(sigy/1e-6)) + "um\n")
        file.write("  RadiusZ = " + str("{:.4f}".format( np.sqrt(2.)*sigz*kp )) + "		# =sqrt(2)sigz*kp - sigx = sigr = " + str("{:.2f}".format(sigz/1e-6)) + "um\n")

        x0 = parameters.get("x01")
        y0 = parameters.get("y01")
        z0 = parameters.get("z01")
        file.write("  x0 = " + str(x0) + "\n")
        file.write("  y0 = " + str(y0) + "			# if 0, on axis injection\n")
        file.write("  z0 = " + str(z0) + "			# if 0, on axis injection\n")
        comment_driver += '(x0,y0,z0)=(' +str(x0)+ ',' + str(y0) + ',' + str(z0) + '); '
        Delta = parameters.get("Delta1")
        file.write("  Delta = " + str(Delta) + "\n")

        E = parameters.get("Energy1")
        comment_driver+='E0=' + str("{:.0f}".format(E)) + 'MeV; '
        Polarity = parameters.get("Polarity1")
        if Polarity == 1:
            gamma = E/E0p+1
        elif Polarity == -1:
            gamma = E/E0+1
        beta = np.sqrt(1.-1./gamma**2)
        Px0 = gamma * beta # longitudinal momentum
        file.write("  Px0 = " + str("{:.4e}".format(Px0)) + "		# normalized momentum - Driver initial energy = "+str("{:.0f}".format(E))+" MeV\n")
        file.write("  Py0 = 0." + "			# no transverse momentum\n")
        file.write("  Pz0 = 0." + "			# no transverse momentum\n")

        P_perCell = parameters.get("P_perCell1")
        comment_driver += 'ppc='+str(int(P_perCell)) + '; '
        file.write("  P_perCell = " + str(int(P_perCell)) + "		# particle per cell: 16, 27, 100\n")

        deltaE = parameters.get("Espread1")
        comment_driver += 'EnergySpread='+str("{:.1f}".format(deltaE*100)) + '%; '
        PspreadX = gamma * beta * deltaE
        file.write("  PspreadX = " + str("{:.4e}".format(PspreadX)) + "		# Depends on energy spread = dE/E = " + str("{:.0f}".format(deltaE*100))+  "%\n")

        epsn = parameters.get("EmittanceN1")
        comment_driver+='en=' + str("{:.1f}".format(epsn/1e-6)) + 'mm-mrad; '
        PspreadY = PspreadZ = epsn/sigr
        file.write("  PspreadY = " + str("{:.4e}".format(PspreadY)) + "		# Depends on emittance (=" + str("{:.2e}".format(epsn)) +" m-rad) and sigr (="+str("{:.2f}".format(sigr/1e-6))+" um)\n")
        file.write("  PspreadZ = " + str("{:.4e}".format(PspreadZ)) + "		# Depends on emittance (=" + str("{:.2e}".format(epsn)) +" m-rad) and sigr (="+str("{:.2f}".format(sigr/1e-6))+" um)\n")

        PhaseSpaceFillFlag = parameters.get("PhaseSpaceFillFlag1")
        file.write("  PhaseSpaceFillFlag = " + str(int(PhaseSpaceFillFlag)) + "\n")
        Polarity = parameters.get("Polarity1")
        file.write("  Polarity = " + str(int(Polarity)) + "			# charge: 1-proton; -1-electron\n")
        if Polarity == -1:
            comment_driver+= " Electron driver; "
        elif Polarity == 1:
            comment_driver+= " Proton driver; "
        MassAE = parameters.get("MassAE1")
        file.write("  MassAE = " + str("{:.4e}".format(MassAE)) + "		# particle mass in atomic units\n")
        Type = parameters.get("Type1")
        file.write("  Type = " + str(int(Type)) + "\n")
        Zombie = parameters.get("Zombie1")
        file.write("  Zombie = " + str(int(Zombie)) + "\n")
        Beam = parameters.get("Beam1")
        file.write("  Beam = " + str(int(Beam)) + "			# if species is beam set to one\n")
        RandomSeed = parameters.get("RandomSeed1")
        file.write("  RandomSeed = " + str(int(RandomSeed)) + "		# 0: grid loading - 1: particles initialised with random position inside cell; when it is beam put it to 1\n")
        SkipSaveFlag = parameters.get("SkipSaveFlag1")
        file.write("  SkipSaveFlag = " + str(int(SkipSaveFlag)) + "		# (=1) to turn off beam particles saving\n")

        file.write("/\n# ************************************************\n")

        if Nspecies == 3: # it means you have plasma, driver and witness
            file.write("&Specie2 			# Witness beam\n")
            #################################################
            Distribution = parameters.get("Distribution2")
            if Distribution == 1:
                comment_witness+= "Gaussian profile; "
            else:
                raise ValueError("check driver distribution. change driver comment.")
            file.write("  Distribution = " + str(int(Distribution)) + "		# Gaussian ellipsoid\n")
            Q = parameters.get("Charge2") * 1e-12 #C
            comment_witness+= 'q='+str("{:.2f}".format(Q/1e-12)) + 'pC; '
            sigx = parameters.get("sigx2")
            comment_witness+= 'sigx='+str("{:.2f}".format(sigx/1e-6)) + 'um; '
            sigr2 = parameters.get("sigr2")
            epsn = parameters.get("EmittanceN2")
            E = parameters.get("Energy2")
            Polarity = parameters.get("Polarity2")
            if Polarity == 1:
                gamma = E/E0p+1
            elif Polarity == -1:
                gamma = E/E0+1
            #gamma = E/E0+1
            beta = np.sqrt(1.-1./gamma**2)
            plasma_matched_radius = ( 2*epsn**2/gamma/kp**2 ) ** (1./4.)
            if sigr2 == 0:
                sigr = plasma_matched_radius
            else:
                sigr = sigr2

            comment_witness+= 'sigr='+str("{:.2f}".format(sigr/1e-6)) + 'um; '
            ne = Q/(2*pi)**(3./2.)/sigx/sigr**2/e
            comment_witness+= 'n/n0='+str("{:.2f}".format(ne/n0)) + '; '
            file.write("  Density = " + str("{:.3e}".format(ne/n0)) + "		# normalized to the plasma density; ne=" + str("{:.2e}".format(ne/1e6)) +"/cm3 and n0=" + str("{:.2e}".format(n0/1e6)) +"/cm3\n")

            Begin = parameters.get("Begin2")
            file.write("  Begin = " + str(Begin) + "\n")
            PlateauBegin = parameters.get("PlateauBegin2")
            file.write("  PlateauBegin = " + str(int(PlateauBegin)) + "\n")
            PlateauEnd = parameters.get("PlateauEnd2")
            file.write("  PlateauEnd = " + str("{:.1e}".format(PlateauEnd)) + "\n")
            End = parameters.get("End2")
            file.write("  End = " + str("{:.1e}".format(End)) + "\n")

            sigy = sigz = sigr
            file.write("  RadiusX = " + str("{:.4f}".format( np.sqrt(2.)*sigx*kp )) + "		# =sqrt(2)sigx*kp - sigx = " + str("{:.2f}".format(sigx/1e-6)) + "um\n")
            file.write("  RadiusY = " + str("{:.4f}".format( np.sqrt(2.)*sigy*kp )) + "		# =sqrt(2)sigy*kp - sigy = sigr = " + str("{:.2f}".format(sigy/1e-6)) + "um (beam matched radius is " + str("{:.2f}".format(plasma_matched_radius/1e-6)) + "um)\n")
            file.write("  RadiusZ = " + str("{:.4f}".format( np.sqrt(2.)*sigz*kp )) + "		# =sqrt(2)sigz*kp - sigx = sigr = " + str("{:.2f}".format(sigz/1e-6)) + "um (beam matched radius is " + str("{:.2f}".format(plasma_matched_radius/1e-6)) + "um)\n")

            x0 = parameters.get("x02")
            y0 = parameters.get("y02")
            z0 = parameters.get("z02")
            comment_witness += '(x0,y0,z0)=(' +str(x0)+ ',' + str(y0) + ',' + str(z0) + '); '
            file.write("  x0 = " + str(x0) + "\n")
            file.write("  y0 = " + str(y0) + "			# if 0, on axis injection\n")
            file.write("  z0 = " + str(z0) + "			# if 0, on axis injection\n")
            Delta = parameters.get("Delta2")
            file.write("  Delta = " + str(Delta) + "\n")

            comment_witness+='E0=' + str("{:.0f}".format(E)) + 'MeV; '
            Px0 = gamma * beta # longitudinal momentum
            file.write("  Px0 = " + str("{:.4e}".format(Px0)) + "		# normalized momentum - Driver initial energy = "+str("{:.0f}".format(E))+" MeV\n")
            file.write("  Py0 = 0." + "			# no transverse momentum\n")
            file.write("  Pz0 = 0." + "			# no transverse momentum\n")

            P_perCell = parameters.get("P_perCell2")
            comment_witness+=str(int(P_perCell)) + '; '
            file.write("  P_perCell = " + str(int(P_perCell)) + "		# particle per cell: 16, 27, 100\n")

            deltaE = parameters.get("Espread2")
            comment_witness += 'EnergySpread='+str("{:.1f}".format(deltaE*100)) + '%; '
            PspreadX = gamma * beta * deltaE
            file.write("  PspreadX = " + str("{:.4e}".format(PspreadX)) + "		# Depends on energy spread = dE/E = " + str("{:.0f}".format(deltaE*100))+  "%\n")

            comment_witness+='en=' + str("{:.1f}".format(epsn/1e-6)) + 'mm-mrad; '
            PspreadY = PspreadZ = epsn/sigr
            file.write("  PspreadY = " + str("{:.4e}".format(PspreadY)) + "		# Depends on emittance (=" + str("{:.2e}".format(epsn)) +" m-rad) and sigr (="+str("{:.2f}".format(sigr/1e-6))+" um)\n")
            file.write("  PspreadZ = " + str("{:.4e}".format(PspreadZ)) + "		# Depends on emittance (=" + str("{:.2e}".format(epsn)) +" m-rad) and sigr (="+str("{:.2f}".format(sigr/1e-6))+" um)\n")

            PhaseSpaceFillFlag = parameters.get("PhaseSpaceFillFlag2")
            file.write("  PhaseSpaceFillFlag = " + str(int(PhaseSpaceFillFlag)) + "\n")
            #Polarity = parameters.get("Polarity2")
            file.write("  Polarity = " + str(int(Polarity)) + "			# charge: 1-proton; -1-electron\n")
            if Polarity == -1:
                comment_witness+= " Electron driver; "
            elif Polarity == 1:
                comment_witness+= " Proton driver; "
            MassAE = parameters.get("MassAE2")
            file.write("  MassAE = " + str("{:.4e}".format(MassAE)) + "		# particle mass in atomic units\n")
            Type = parameters.get("Type2")
            file.write("  Type = " + str(int(Type)) + "\n")
            Zombie = parameters.get("Zombie2")
            file.write("  Zombie = " + str(int(Zombie)) + "\n")
            Beam = parameters.get("Beam2")
            file.write("  Beam = " + str(int(Beam)) + "			# if species is beam set to one\n")
            RandomSeed = parameters.get("RandomSeed2")
            file.write("  RandomSeed = " + str(int(RandomSeed)) + "		# 0: grid loading - 1: particles initialised with random position inside cell; when it is beam put it to 1\n")
            SkipSaveFlag = parameters.get("SkipSaveFlag2")
            file.write("  SkipSaveFlag = " + str(int(SkipSaveFlag)) + "		# (=1) to turn off beam particles saving\n")
            file.write("/\n# ************************************************\n")

        file.write("&Synchrotron			# Scan photons within the energy Emin-Emax\n")
        ############################
        Emax = parameters.get("Emax")
        file.write("  Emax = " + str("{:.1e}".format(Emax)) + "		# [eV]\n")

        Emin = parameters.get("Emin")
        file.write("  Emin = " + str("{:.1e}".format(Emin)) + "		# [eV]\n")

        comment_radiation+='Emin=' + str("{:.1f}".format(Emin)) + 'eV; '
        comment_radiation+='Emax=' + str("{:.1f}".format(Emax/1000)) + 'keV; '

        ThetaMax = parameters.get("ThetaMax")
        comment_radiation+='ThetaMax=' + str("{:.1f}".format(ThetaMax*1000)) + 'mrad; '
        file.write("  ThetaMax = " + str("{:.3f}".format(ThetaMax)) + "		# rad; max theta\n")

        nPhibins = parameters.get("nPhibins")
        file.write("  nPhibins = " + str(int(nPhibins)) + "		# number of bins to divide phi angle\n")

        nThetabins = parameters.get("nThetabins")
        file.write("  nThetabins = " + str(int(nThetabins)) + "		# number of bins to divide theta angle\n")

        nEbins = parameters.get("nEbins")
        file.write("  nEbins = " + str(int(nEbins)) + "			# number of bins to divide E\n")

        comment_radiation += '(nE,nTheta,nPhi)=(' + str(int(nEbins)) + ',' + str(int(nThetabins)) + ',' + str(int(nPhibins)) + '); '
        SynMin = parameters.get("SynMin")
        comment_radiation+= 'minPhotonEnergy=' + str("{:.1f}".format(SynMin)) + 'eV; '
        file.write("  SynMin = " + str("{:.1e}".format(SynMin)) + "		# [eV] minimum photon energy to be counted\n")

        phase = PhaseStop_rounded
        file.write("  phase = " + str("{:.1f}".format(phase)) + "		# (=PhaseStop=" + str("{:.0f}".format(PhaseStop_rounded)) +") maximum phase of calculation;  \n")

        screen_to_cell_distance = parameters.get("screen_to_cell_distance")
        ScreenPosition = PhaseStop_rounded + screen_to_cell_distance*kp
        comment_radiation = '# Screen is ' + str(screen_to_cell_distance) + ' m after plasma cell.'
        file.write("  ScreenPosition = " + str("{:.1f}".format(ScreenPosition)) +  "      " +comment_radiation + "\n")

        file.write("/\n# ************************************************\n")

        file.write("&Movie2dHDF5			# Parameters to save in 2D\n")
        ############################
        Frame0 = parameters.get("Frame0")
        Frame1 = parameters.get("Frame1")
        Frame2 = parameters.get("Frame2")
        Frame3 = parameters.get("Frame3")
        Frame4 = parameters.get("Frame4")
        Frame5 = parameters.get("Frame5")

        file.write("  Frame0 = " + Frame0 + "\n")
        file.write("  Frame1 = " + str("n1") + "\n")
        file.write("  Frame2 = " + str("n2") + "\n")
        file.write("  Frame3 = " + str("ex") + "\n")
        file.write("  Frame4 = " + str("ForceY") + "		# ForcY = Ey-c*Bz\n")
        file.write("  Frame5 = " + str("psi") + "\n")

        file.write("/\n# ************************************************\n")

        file.write("&Boundary_Xm\n")
        FieldCondition = parameters.get("FieldCondition_Boundary_Xm")
        ParticlesCondition = parameters.get("ParticlesCondition_Boundary_Xm")
        file.write("  FieldCondition = "+ str(FieldCondition) +"\n")
        file.write("  ParticlesCondition = "+ str(ParticlesCondition) +"\n")
        file.write("/\n")


        file.write("&Boundary_Xp\n")
        FieldCondition = parameters.get("FieldCondition_Boundary_Xp")
        ParticlesCondition = parameters.get("ParticlesCondition_Boundary_Xp")
        file.write("  FieldCondition = "+ str(FieldCondition) +"\n")
        file.write("  ParticlesCondition = "+ str(ParticlesCondition) +"\n")
        file.write("/\n")

        file.write("&Boundary_Ym\n")
        FieldCondition = parameters.get("FieldCondition_Boundary_Ym")
        ParticlesCondition = parameters.get("ParticlesCondition_Boundary_Ym")
        file.write("  FieldCondition = "+ str(FieldCondition) +"\n")
        file.write("  ParticlesCondition = "+ str(ParticlesCondition) +"\n")
        file.write("/\n")

        file.write("&Boundary_Yp\n")
        FieldCondition = parameters.get("FieldCondition_Boundary_Yp")
        ParticlesCondition = parameters.get("ParticlesCondition_Boundary_Yp")
        file.write("  FieldCondition = "+ str(FieldCondition) +"\n")
        file.write("  ParticlesCondition = "+ str(ParticlesCondition) +"\n")
        file.write("/\n")

        file.write("&Boundary_Zm\n")
        FieldCondition = parameters.get("FieldCondition_Boundary_Zm")
        ParticlesCondition = parameters.get("ParticlesCondition_Boundary_Zm")
        file.write("  FieldCondition = "+ str(FieldCondition) +"\n")
        file.write("  ParticlesCondition = "+ str(ParticlesCondition) +"\n")
        file.write("/\n")

        file.write("&Boundary_Zp\n")
        FieldCondition = parameters.get("FieldCondition_Boundary_Zp")
        ParticlesCondition = parameters.get("ParticlesCondition_Boundary_Zp")
        file.write("  FieldCondition = "+ str(FieldCondition) +"\n")
        file.write("  ParticlesCondition = "+ str(ParticlesCondition) +"\n")
        file.write("/\n")
        file.close()

    # Open the file for reading and writing
    with open(input_deck, "r+") as file:
    # Read the contents of the file
        lines = file.readlines()

        # Insert the new line after the 5th line
        lines.insert(5, comment_plasma+"\n")
        lines.insert(6, comment_driver+"\n")
        lines.insert(7, comment_witness+"\n")
        lines.insert(8, comment_simulation+"\n")
        lines.insert(9, comment_radiation+"\n")
        lines.insert(10, '# ************************************************\n')


        # Move the file cursor to the beginning of the file
        file.seek(0)

        # Write the modified contents back to the file
        file.writelines(lines)
        file.close()


def main(input_file):
    parameters = read_default_parameters(input_file)
    generat_inputdeck(parameters, 'v.ini')

if __name__ == "__main__":
    # Check if three arguments are provided
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_file>")
    else:
        # Parse command-line arguments
        input_file = sys.argv[1]

        # Call the main function with the arguments
        main(input_file)
