/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.prediction.homology.gui;


import javax.swing.*;
import java.awt.*;

public class InflateParametersDialog extends JDialog {
    public JLabel labelInflateEnergyWeight;
    public JTextField textInflateEnergyWeight;
    public JLabel labelInflateEnergyRmsTarget;
    public JTextField textInflateEnergyRmsTarget;
    public JButton doneInflateParameters;

    public InflateParametersDialog() {
        InflateParametersLayout customLayout = new InflateParametersLayout();

        getContentPane().setFont(new Font("Helvetica", Font.PLAIN, 12));
        getContentPane().setLayout(customLayout);

        labelInflateEnergyWeight = new JLabel("inflateEnergy weight ");
        getContentPane().add(labelInflateEnergyWeight);

        textInflateEnergyWeight = new JTextField("1.2");
        getContentPane().add(textInflateEnergyWeight);

        labelInflateEnergyRmsTarget = new JLabel("inflateEnergy RmsTarget ");
        getContentPane().add(labelInflateEnergyRmsTarget);

        textInflateEnergyRmsTarget = new JTextField("1.5");
        getContentPane().add(textInflateEnergyRmsTarget);

        doneInflateParameters = new JButton("Done");
        getContentPane().add(doneInflateParameters);

        pack();


    }


}

class InflateParametersLayout implements LayoutManager {

    public InflateParametersLayout() {
    }

    public void addLayoutComponent(String name, Component comp) {
    }

    public void removeLayoutComponent(Component comp) {
    }

    public Dimension preferredLayoutSize(Container parent) {
        Dimension dim = new Dimension(0, 0);

        Insets insets = parent.getInsets();
        dim.width = 331 + insets.left + insets.right;
        dim.height = 138 + insets.top + insets.bottom;

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
            c.setBounds(insets.left + 8, insets.top + 8, 232, 24);
        }
        c = parent.getComponent(1);
        if (c.isVisible()) {
            c.setBounds(insets.left + 248, insets.top + 8, 72, 24);
        }
        c = parent.getComponent(2);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 32, 232, 24);
        }
        c = parent.getComponent(3);
        if (c.isVisible()) {
            c.setBounds(insets.left + 248, insets.top + 32, 72, 24);
        }
        c = parent.getComponent(4);
        if (c.isVisible()) {
            c.setBounds(insets.left + 208, insets.top + 96, 72, 24);
        }
    }
}
