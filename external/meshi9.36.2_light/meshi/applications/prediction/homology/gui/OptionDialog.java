/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.prediction.homology.gui;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public class OptionDialog extends JDialog {
    public final JLabel parametersExplanationLabel;
    public final JLabel numberOfModelsLabel;
    public final JLabel labelSeedNumber;
    public final JTextField numberOfModelsText;
    public final JButton buttonDistanceConstrainParameters;
    public final JButton buttonMinimizationParametersShort;
    public final JButton buttonMinimizationParametersFull;
    public final JButton buttonDoneOptionDialog;
    public final JTextField textSeedNumber;
    public final JButton buttonEnergyWeights;
    public final JButton buttonInflate;
    public final JButton buttonParametersForChainCompletion;
    public final JButton buttonNonFrozenRange;

    public OptionDialog() {
        OptionDialog3Layout customLayout = new OptionDialog3Layout();

        getContentPane().setFont(new Font("Helvetica", Font.PLAIN, 12));
        getContentPane().setLayout(customLayout);

        parametersExplanationLabel = new JLabel("Choose your prefered commands and parameters");
        getContentPane().add(parametersExplanationLabel);

        numberOfModelsLabel = new JLabel("Number Of Models");
        getContentPane().add(numberOfModelsLabel);

        labelSeedNumber = new JLabel("Seed Number");
        getContentPane().add(labelSeedNumber);

        numberOfModelsText = new JTextField("1");
        getContentPane().add(numberOfModelsText);

        buttonDistanceConstrainParameters = new JButton("Distance Constrain Parameters");
        getContentPane().add(buttonDistanceConstrainParameters);

        buttonMinimizationParametersShort = new JButton("Minimization parameters - short relaxation");
        getContentPane().add(buttonMinimizationParametersShort);

        buttonMinimizationParametersFull = new JButton("Minimization parameters - full minimization");
        getContentPane().add(buttonMinimizationParametersFull);

        buttonDoneOptionDialog = new JButton("Done");
        getContentPane().add(buttonDoneOptionDialog);

        textSeedNumber = new JTextField("1");
        getContentPane().add(textSeedNumber);

        buttonEnergyWeights = new JButton("Energy Weights");
        getContentPane().add(buttonEnergyWeights);

        buttonInflate = new JButton("InflateBySegments parameters");
        getContentPane().add(buttonInflate);

        buttonParametersForChainCompletion = new JButton("Parameters for Chain Completion");
        getContentPane().add(buttonParametersForChainCompletion);

        buttonNonFrozenRange = new JButton("NonFrozen Range");
        getContentPane().add(buttonNonFrozenRange);

        setSize(getPreferredSize());

        addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                System.exit(0);
            }
        });
    }

    public static void main(String args[]) {
        OptionDialog window = new OptionDialog();

        window.setTitle("OptionDialog3");
        window.pack();
        window.show();
    }
}

class OptionDialog3Layout implements LayoutManager {

    public OptionDialog3Layout() {
    }

    public void addLayoutComponent(String name, Component comp) {
    }

    public void removeLayoutComponent(Component comp) {
    }

    public Dimension preferredLayoutSize(Container parent) {
        Dimension dim = new Dimension(0, 0);

        Insets insets = parent.getInsets();
        dim.width = 500 + insets.left + insets.right;
        dim.height = 334 + insets.top + insets.bottom;

        return dim;
    }

    public Dimension minimumLayoutSize(Container parent) {
        Dimension dim = new Dimension(0, 0);
        return dim;
    }

    public void layoutContainer(Container parent) {
        Insets insets = parent.getInsets();

        Component c;
        c = parent.getComponent(0);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 8, 480, 32);
        }
        c = parent.getComponent(1);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 56, 152, 24);
        }
        c = parent.getComponent(2);
        if (c.isVisible()) {
            c.setBounds(insets.left + 256, insets.top + 56, 152, 24);
        }
        c = parent.getComponent(3);
        if (c.isVisible()) {
            c.setBounds(insets.left + 160, insets.top + 56, 72, 24);
        }
        c = parent.getComponent(4);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 200, 408, 24);
        }
        c = parent.getComponent(5);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 176, 408, 24);
        }
        c = parent.getComponent(6);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 152, 408, 24);
        }
        c = parent.getComponent(7);
        if (c.isVisible()) {
            c.setBounds(insets.left + 264, insets.top + 288, 72, 24);
        }
        c = parent.getComponent(8);
        if (c.isVisible()) {
            c.setBounds(insets.left + 408, insets.top + 56, 24, 24);
        }
        c = parent.getComponent(9);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 104, 408, 24);
        }
        c = parent.getComponent(10);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 128, 408, 24);
        }
        c = parent.getComponent(11);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 224, 408, 24);
        }
        c = parent.getComponent(12);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 248, 408, 24);
        }
    }
}
