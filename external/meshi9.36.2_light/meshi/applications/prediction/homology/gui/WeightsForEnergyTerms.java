/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.prediction.homology.gui;

import java.awt.*;
import javax.swing.*;

public class WeightsForEnergyTerms extends JFrame {
    public JLabel labelBondedEnergyWeight;
    public JTextField textBondedEnergyWeight;
    public JLabel labelTwoTorsionsEnergyWeight;
    public JLabel labelPlaneEnergyWeight;
    public JLabel labelAngleEnergyWeight;
    public JLabel labelLennardJonesWeight;
    public JLabel labelLennardJonesCaWeight;
    public JTextField textAngleEnergyWeight;
    public JTextField textLennardJonesWeight;
    public JTextField textLennardJonesCaWeight;
    public JTextField textTwoTorsionsEnergyWeight;
    public JTextField textPlaneEnergyWeight;
    public JLabel labelSolvateCA100EnergyWeight;
    public JLabel labelPropensityTorsionEnergyWeight;
    public JLabel labelFlatRamachEnergyWeight;
    public JTextField textSolvateCA100EnergyWeight;
    public JTextField textPropensityTorsionEnergyWeight;
    public JTextField textFlatRamachEnergyWeight;
    public JButton doneWeightsForEnergyTerms;
    public JLabel labelSolventBb60;
    public JTextField textSolventBb60;
    public JTextField textSolventAa40;
    public JLabel labelSolventAa40;
    public JLabel labelOutOfPlaneEnergyWeight;
    public JTextField textOutOfPlaneEnergyWeight;
    public JLabel labelAlphaAngleEnergyWeight;
    public JTextField textAlphaAngleEnergyWeight;
    public JLabel labelAlphaTorsionEnergyWeight;
    public JTextField textAlphaTorsionEnergyWeight;
    public JLabel labelHydrogenBondsWeight;
    public JTextField textHydrogenBondsWeight;
    public JLabel labelHydrogenBondsPairsWeight;
    public JTextField textHydrogenBondsPairsWeight;
    public JLabel labelExcludedVolumeWeight;
    public JTextField textExcludedVolumeWeight;

    public WeightsForEnergyTerms() {
        WeightsForEnergyTermsLayout customLayout = new WeightsForEnergyTermsLayout();

        getContentPane().setFont(new Font("Helvetica", Font.PLAIN, 12));
        getContentPane().setLayout(customLayout);

        labelBondedEnergyWeight = new JLabel("bondEnergy weight");
        getContentPane().add(labelBondedEnergyWeight);

        textBondedEnergyWeight = new JTextField("1.0");
        getContentPane().add(textBondedEnergyWeight);

        labelTwoTorsionsEnergyWeight = new JLabel("twoTorsionsEnergy weight");
        getContentPane().add(labelTwoTorsionsEnergyWeight);

        labelPlaneEnergyWeight = new JLabel("planeEnergy weight");
        getContentPane().add(labelPlaneEnergyWeight);

        labelAngleEnergyWeight = new JLabel("angleEnergy weight");
        getContentPane().add(labelAngleEnergyWeight);

        labelLennardJonesWeight = new JLabel("unsupportedEnergies weight");
        getContentPane().add(labelLennardJonesWeight);

        labelLennardJonesCaWeight = new JLabel("LennardJonesCa weight");
        getContentPane().add(labelLennardJonesCaWeight);

        textAngleEnergyWeight = new JTextField("1.0");
        getContentPane().add(textAngleEnergyWeight);

        textLennardJonesWeight = new JTextField("0.2");
        getContentPane().add(textLennardJonesWeight);

        textLennardJonesCaWeight = new JTextField("0.2");
        getContentPane().add(textLennardJonesCaWeight);

        textTwoTorsionsEnergyWeight = new JTextField("1.0");
        getContentPane().add(textTwoTorsionsEnergyWeight);

        textPlaneEnergyWeight = new JTextField("1.0");
        getContentPane().add(textPlaneEnergyWeight);

        labelSolvateCA100EnergyWeight = new JLabel("solvateCA100Energy weight");
        getContentPane().add(labelSolvateCA100EnergyWeight);

        labelPropensityTorsionEnergyWeight = new JLabel("propensityTorsionEnergy weight");
        getContentPane().add(labelPropensityTorsionEnergyWeight);

        labelFlatRamachEnergyWeight = new JLabel("flatRamachEnergy weight");
        getContentPane().add(labelFlatRamachEnergyWeight);

        textSolvateCA100EnergyWeight = new JTextField("1.0");
        getContentPane().add(textSolvateCA100EnergyWeight);

        textPropensityTorsionEnergyWeight = new JTextField("1.0");
        getContentPane().add(textPropensityTorsionEnergyWeight);

        textFlatRamachEnergyWeight = new JTextField("10.0");
        getContentPane().add(textFlatRamachEnergyWeight);

        doneWeightsForEnergyTerms = new JButton("Done");
        getContentPane().add(doneWeightsForEnergyTerms);

        labelSolventBb60 = new JLabel("solvateBB60Energy weight ");
        getContentPane().add(labelSolventBb60);

        textSolventBb60 = new JTextField("1.0");
        getContentPane().add(textSolventBb60);

        textSolventAa40 = new JTextField("1.0");
        getContentPane().add(textSolventAa40);

        labelSolventAa40 = new JLabel("solvateAA40Energy weight ");
        getContentPane().add(labelSolventAa40);

        labelOutOfPlaneEnergyWeight = new JLabel("outOfPlaneEnergy weight ");
        getContentPane().add(labelOutOfPlaneEnergyWeight);

        textOutOfPlaneEnergyWeight = new JTextField("1.0");
        getContentPane().add(textOutOfPlaneEnergyWeight);

        labelAlphaAngleEnergyWeight = new JLabel("alphaAngleEnergy weight");
        getContentPane().add(labelAlphaAngleEnergyWeight);

        textAlphaAngleEnergyWeight = new JTextField("1.0");
        getContentPane().add(textAlphaAngleEnergyWeight);

        labelAlphaTorsionEnergyWeight = new JLabel("alphaTorsionEnergy weight ");
        getContentPane().add(labelAlphaTorsionEnergyWeight);

        textAlphaTorsionEnergyWeight = new JTextField("1.0");
        getContentPane().add(textAlphaTorsionEnergyWeight);

        labelHydrogenBondsWeight = new JLabel("hydrogenBonds weight ");
        getContentPane().add(labelHydrogenBondsWeight);

        textHydrogenBondsWeight = new JTextField("0.0");
        getContentPane().add(textHydrogenBondsWeight);

        labelHydrogenBondsPairsWeight = new JLabel("hydrogenBondsPairs weight ");
        getContentPane().add(labelHydrogenBondsPairsWeight);

        textHydrogenBondsPairsWeight = new JTextField("0.0");
        getContentPane().add(textHydrogenBondsPairsWeight);

        labelExcludedVolumeWeight = new JLabel("excludedVolume weight ");
        getContentPane().add(labelExcludedVolumeWeight);

        textExcludedVolumeWeight = new JTextField("1.0");
        getContentPane().add(textExcludedVolumeWeight);

        pack();

    }


}

