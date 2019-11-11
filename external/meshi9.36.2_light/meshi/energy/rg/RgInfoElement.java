/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.rg;

import meshi.dataStructures.Array;
import meshi.energy.EnergyInfoElement;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfo;

import java.util.ArrayList;

public class RgInfoElement {//extends EnergyInfoElement {


    public final double all, RGlogaritmicWeight, RGratioWeight, RGpolarWeight, RGnonPolarWeight, RGbackboneWeight;
    public final EnergyInfoElement N_RG_HSS, E_HSS_HCOIL, E_HSS_BSS, E_HSS_BCOIL, E_HSS_PSS, E_HSS_PCOIL, HSS, BSS, PSS, HSS_HCOIL, HSS_BSS, HSS_BCOIL, HSS_PSS, HSS_PCOIL;

    public RgInfoElement(double all, double RGlogaritmicWeight, double RGratioWeight, double RGpolarWeight, double RGnonPolarWeight, double RGbackboneWeight) {
        //super(InfoType.RG, "Radius of Gyration related terms ");
        this.all = all;
        this.RGlogaritmicWeight = RGlogaritmicWeight;
        this.RGratioWeight = RGratioWeight;
        this.RGpolarWeight = RGpolarWeight;
        this.RGnonPolarWeight = RGnonPolarWeight;
        this.RGbackboneWeight = RGbackboneWeight;
        ArrayList<MeshiInfo> infoList = new ArrayList();
        infoList.add(N_RG_HSS    = new EnergyInfoElement(InfoType.N_RG_HSS, "log(Length)  Vs. log(RG of hydrophobic side-chains in secondary structure elements", Double.MIN_VALUE));
        infoList.add(E_HSS_HCOIL = new EnergyInfoElement(InfoType.E_HSS_HCOIL, "log(RG of hydrophobic side-chains in secondary structure elements Vs of hydrophobic side-chains in coil", Double.MIN_VALUE));
        infoList.add(E_HSS_BSS   = new EnergyInfoElement(InfoType.E_HSS_BSS, "log(RG of hydrophobic side-chains in secondary structure elements Vs CAs in SSEs", Double.MIN_VALUE));
        infoList.add(E_HSS_BCOIL = new EnergyInfoElement(InfoType.E_HSS_BCOIL, "log(RG of hydrophobic side-chains in secondary structure elements Vs CAs in coil", Double.MIN_VALUE));
        infoList.add(E_HSS_PSS   = new EnergyInfoElement(InfoType.E_HSS_CSS, "log(RG of hydrophobic side-chains in secondary structure elements Vs Polar side chain getAtoms in SSEs", Double.MIN_VALUE));
        infoList.add(E_HSS_PCOIL = new EnergyInfoElement(InfoType.E_HSS_CCOIL, "log(RG of hydrophobic side-chains in secondary structure elements Vs Polar side chain getAtoms in Coil", Double.MIN_VALUE));
        infoList.add(HSS         = new EnergyInfoElement(InfoType.HSS, "log(RG of hydrophobic side-chains in secondary structure elements", Double.MIN_VALUE));
        infoList.add(BSS         = new EnergyInfoElement(InfoType.BSS, "log(RG of CAs in secondary structure elements", Double.MIN_VALUE));
        infoList.add(PSS         = new EnergyInfoElement(InfoType.CSS, "log(RG of polar getAtoms in secondary structure elements", Double.MIN_VALUE));
        infoList.add(HSS_HCOIL   = new EnergyInfoElement(InfoType.HSS_HCOIL, "log(RG of hydrophobic side-chains in secondary structure elements Vs of hydrophobic side-chains in coil", Double.MIN_VALUE));
        infoList.add(HSS_BSS     = new EnergyInfoElement(InfoType.HSS_BSS, "log(RG of hydrophobic side-chains in secondary structure elements Vs CAs in SSEs", Double.MIN_VALUE));
        infoList.add(HSS_BCOIL   = new EnergyInfoElement(InfoType.HSS_BCOIL, "log(RG of hydrophobic side-chains in secondary structure elements Vs CAs in coil", Double.MIN_VALUE));
        infoList.add(HSS_PSS     = new EnergyInfoElement(InfoType.HSS_CSS, "log(RG of hydrophobic side-chains in secondary structure elements Vs Polar side chain getAtoms in SSEs", Double.MIN_VALUE));
        infoList.add(HSS_PCOIL   = new EnergyInfoElement(InfoType.HSS_CCOIL, "log(RG of hydrophobic side-chains in secondary structure elements Vs Polar side chain getAtoms in Coil", Double.MIN_VALUE));
    }
}