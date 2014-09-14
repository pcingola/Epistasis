package meshi.energy.solvate;
import meshi.molecularElements.atoms.*;
import meshi.util.*;
import meshi.geometry.*;
import meshi.parameters.*;
import meshi.molecularElements.hydrogenBonds.*;
import meshi.util.mathTools.*;

public class SolvateDistanceAttributeCreator {
    public static SolvateDistanceAttribute create(Distance dis, SolvateParametersList parameters, AbstractHydrogenBondList solvateHB) {
	 // Does this distance involve a hydrogen ??
	 // ----------------------------------------
	boolean isDisInvolvesHydrogens;
	 if (dis.atom1().type().isHydrogen() || dis.atom2().type().isHydrogen()) 
	     return new SolvateDistanceAttribute(true, false);
	 else
	     isDisInvolvesHydrogens = false;

	 // Does this distance involves two polar atoms ??
	 // ----------------------------------------------
	 boolean isDisBetween2Polars;
	 if ((dis.atom1().type().isOxygen() || dis.atom1().type().isNitrogen() || (dis.atom1().type() == AtomType.CSG)) && 
	     (dis.atom2().type().isOxygen() || dis.atom2().type().isNitrogen() || (dis.atom2().type() == AtomType.CSG))) 
	     isDisBetween2Polars = true;
	 else 
	     isDisBetween2Polars = false;

	 //---------------------
	double saltBridgeFactorForHBenergyA1, saltBridgeFactorForHBenergyA2, saltBridgeFactorA1, saltBridgeFactorA2;
    	if  ((((dis.atom1().type() == AtomType.KNZ) || (dis.atom1().type() == AtomType.RNH) || (dis.atom1().type() == AtomType.TRN)) && 
	      ((dis.atom2().type() == AtomType.DOD) || (dis.atom2().type() == AtomType.EOE) || (dis.atom2().type() == AtomType.TRO)))      || 
    	     (((dis.atom2().type() == AtomType.KNZ) || (dis.atom2().type() == AtomType.RNH) || (dis.atom2().type() == AtomType.TRN)) && 
	      ((dis.atom1().type() == AtomType.DOD) || (dis.atom1().type() == AtomType.EOE) || (dis.atom1().type() == AtomType.TRO)))) {  // This is a salt bridge
	    saltBridgeFactorForHBenergyA1 = SolvateEnergy.SALT_BRIDGE_STRENGTH_GENERAL;
	    saltBridgeFactorForHBenergyA2 = SolvateEnergy.SALT_BRIDGE_STRENGTH_GENERAL;
	    if (dis.atom1().type() == AtomType.DOD)
    	     	saltBridgeFactorA1 = SolvateEnergy.SALT_BRIDGE_STRENGTH_ASP_OD;
	    else 
		if (dis.atom1().type() == AtomType.EOE)
		    saltBridgeFactorA1 = SolvateEnergy.SALT_BRIDGE_STRENGTH_GLU_OE;
		else 
		    if (dis.atom1().type() == AtomType.TRO)
			saltBridgeFactorA1 = SolvateEnergy.SALT_BRIDGE_STRENGTH_TRO;
		    else 
			if (dis.atom1().type() == AtomType.KNZ)
			    saltBridgeFactorA1 = SolvateEnergy.SALT_BRIDGE_STRENGTH_LYS_NZ;
			else 
			    if (dis.atom1().type() == AtomType.RNH)
				saltBridgeFactorA1 = SolvateEnergy.SALT_BRIDGE_STRENGTH_ARG_NH;    		     
			    else 
				if (dis.atom1().type() == AtomType.TRN)
				    saltBridgeFactorA1 = SolvateEnergy.SALT_BRIDGE_STRENGTH_TRN; 
				else saltBridgeFactorA1 = 1.0;
	    if (dis.atom2().type() == AtomType.DOD)
		saltBridgeFactorA2 = SolvateEnergy.SALT_BRIDGE_STRENGTH_ASP_OD;
	    else 
		if (dis.atom2().type() == AtomType.EOE)
		    saltBridgeFactorA2 = SolvateEnergy.SALT_BRIDGE_STRENGTH_GLU_OE;
		else 
		    if (dis.atom2().type() == AtomType.TRO)
			saltBridgeFactorA2 = SolvateEnergy.SALT_BRIDGE_STRENGTH_TRO;
		    else 
			if (dis.atom2().type() == AtomType.KNZ)
			    saltBridgeFactorA2 = SolvateEnergy.SALT_BRIDGE_STRENGTH_LYS_NZ;
			else 
			    if (dis.atom2().type() == AtomType.RNH)
				saltBridgeFactorA2 = SolvateEnergy.SALT_BRIDGE_STRENGTH_ARG_NH;    		     
			    else 
				if (dis.atom2().type() == AtomType.TRN)
				    saltBridgeFactorA2 = SolvateEnergy.SALT_BRIDGE_STRENGTH_TRN;    
				else saltBridgeFactorA2 = 1.0;
    	}
	else {
	    saltBridgeFactorForHBenergyA1 = 1.0;
	    saltBridgeFactorForHBenergyA2 = 1.0;
	    saltBridgeFactorA1 = 1.0;
	    saltBridgeFactorA2 = 1.0;
	}
	if (isDisBetween2Polars) return new SolvateDistanceAttributeBetweenPolars(dis, solvateHB, 
										  saltBridgeFactorForHBenergyA1, saltBridgeFactorForHBenergyA2, 
										  saltBridgeFactorA1, saltBridgeFactorA2);
	 //------------------------------------
	int TsaiAtomicType1 = parameters.atomicTypeConverter[dis.atom1().type().ordinal()]; // Converting from the 190 atom types to the 14 defined in Tsai 99'
	int TsaiAtomicType2 = parameters.atomicTypeConverter[dis.atom2().type().ordinal()]; // Converting from the 190 atom types to the 14 defined in Tsai 99'
	//=================================
	// Is atom2 non-polar? If so it should contribute to atom1's CNC
	double Cend_12, Cp1_12, Cp2_12, CvalAtp1_12, CvalAtp2_12;
    	if (dis.atom2.type().isCarbon() || (dis.atom2.type() == AtomType.MSD)) {
	    Cend_12     = parameters.Cend[TsaiAtomicType1][TsaiAtomicType2];
	    Cp1_12      = parameters.Cp1[TsaiAtomicType1][TsaiAtomicType2];
	    Cp2_12      = parameters.Cp2[TsaiAtomicType1][TsaiAtomicType2];
	    CvalAtp1_12 = parameters.CvalAtp1[TsaiAtomicType1][TsaiAtomicType2];
	    CvalAtp2_12 = parameters.CvalAtp2[TsaiAtomicType1][TsaiAtomicType2];
	}
	else Cend_12 = Cp1_12 = Cp2_12 = CvalAtp1_12 = CvalAtp2_12 = 0;
       // Is atom1 non-polar? If so it should contribute to atom2's CNC
	double Cend_21, Cp1_21, Cp2_21, CvalAtp1_21, CvalAtp2_21;
        if (dis.atom1.type().isCarbon() || (dis.atom1.type() == AtomType.MSD)) {
            Cend_21     = parameters.Cend[TsaiAtomicType2][TsaiAtomicType1];
            Cp1_21      = parameters.Cp1[TsaiAtomicType2][TsaiAtomicType1];
            Cp2_21      = parameters.Cp2[TsaiAtomicType2][TsaiAtomicType1];
            CvalAtp1_21 = parameters.CvalAtp1[TsaiAtomicType2][TsaiAtomicType1];
            CvalAtp2_21 = parameters.CvalAtp2[TsaiAtomicType2][TsaiAtomicType1];
        }
        else Cend_21 = Cp1_21 = Cp2_21 = CvalAtp1_21 = CvalAtp2_21 = 0;
 
	 return new SolvateDistanceAttributeWithNonPolar(dis, dis.atom1, dis.atom2,  dis.atom1.number, dis.atom2.number, 
					     isDisInvolvesHydrogens, isDisBetween2Polars,
					     TsaiAtomicType1, TsaiAtomicType2,
					     Cend_12, Cp1_12, Cp2_12, CvalAtp1_12, CvalAtp2_12,
					     Cend_21, Cp1_21, Cp2_21, CvalAtp1_21, CvalAtp2_21);
    }
}