/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.prediction.homology.gui;


import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public class HomologyFrame extends JFrame {
    public final JLabel meshiExplanationLabel;
    public final JTextField targetFilePathText;
    public final JTextField alnFiltPathText;
    public final JTextField outputFilePathText;
    public final JButton submit;
    public final JButton commandsButton;
    public final JLabel statusLabel;
    public final JButton targetFileButton;
    public final JButton alnFileButton;
    public final JTextField templateFilePathText;
    public final JButton templateFileButton;
    public final JButton outputDirectoryButton;
    public final JTextField textSecondaryStructure;
    public final JButton buttonSecondaryStructure;
    public final JTextField textTargetName;
    public final JLabel labelTargetName;
    public final JButton buttonReset;
    public final JButton buttonExitSave;
    public final JButton buttonParameters;
    public final JTextField textfieldParameters;
    public final JLabel labelParameters;

    public HomologyFrame() {
        HomologyFrame11Layout customLayout = new HomologyFrame11Layout();

        getContentPane().setFont(new Font("Helvetica", Font.PLAIN, 12));
        getContentPane().setLayout(customLayout);

        meshiExplanationLabel = new JLabel("Please insert your files path and Options then click submit button");
        getContentPane().add(meshiExplanationLabel);

        targetFilePathText = new JTextField("E:\\1.17\\meshi.1.17.export\\examples\\exTarget.txt");
        getContentPane().add(targetFilePathText);

        alnFiltPathText = new JTextField("E:\\1.17\\meshi.1.17.export\\examples\\exAln.txt");
        getContentPane().add(alnFiltPathText);

        outputFilePathText = new JTextField("");
        getContentPane().add(outputFilePathText);

        submit = new JButton("submit");
        getContentPane().add(submit);

        commandsButton = new JButton("Options");
        getContentPane().add(commandsButton);

        statusLabel = new JLabel("");
        getContentPane().add(statusLabel);

        targetFileButton = new JButton("Target");
        getContentPane().add(targetFileButton);

        alnFileButton = new JButton("Alignment");
        getContentPane().add(alnFileButton);

        templateFilePathText = new JTextField("E:\\1.17\\meshi.1.17.export\\examples\\ex.pdb");
        getContentPane().add(templateFilePathText);

        templateFileButton = new JButton("Template");
        getContentPane().add(templateFileButton);

        outputDirectoryButton = new JButton("Output Directory");
        getContentPane().add(outputDirectoryButton);

        textSecondaryStructure = new JTextField("");
        getContentPane().add(textSecondaryStructure);

        buttonSecondaryStructure = new JButton("Secondary Structure");
        getContentPane().add(buttonSecondaryStructure);

        textTargetName = new JTextField("target");
        getContentPane().add(textTargetName);

        labelTargetName = new JLabel("Enter your target protein name");
        getContentPane().add(labelTargetName);

        buttonReset = new JButton("Reset Options");
        getContentPane().add(buttonReset);

        buttonExitSave = new JButton("Exit");
        getContentPane().add(buttonExitSave);

        buttonParameters = new JButton("Parameters Directory");
        getContentPane().add(buttonParameters);

        textfieldParameters = new JTextField("");
        getContentPane().add(textfieldParameters);

        labelParameters = new JLabel("please insert the location for the parameters directory");
        getContentPane().add(labelParameters);

        setSize(getPreferredSize());

        addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                System.exit(0);
            }
        });
    }

    public static void main(String args[]) {
        HomologyFrame window = new HomologyFrame();

        window.setTitle("HomologyFrame11");
        window.pack();
        window.show();
    }
}

class HomologyFrame11Layout implements LayoutManager {

    public HomologyFrame11Layout() {
    }

    public void addLayoutComponent(String name, Component comp) {
    }

    public void removeLayoutComponent(Component comp) {
    }

    public Dimension preferredLayoutSize(Container parent) {
        Dimension dim = new Dimension(0, 0);

        Insets insets = parent.getInsets();
        dim.width = 652 + insets.left + insets.right;
        dim.height = 517 + insets.top + insets.bottom;

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
            c.setBounds(insets.left + 8, insets.top + 120, 608, 16);
        }
        c = parent.getComponent(1);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 152, 416, 24);
        }
        c = parent.getComponent(2);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 192, 416, 24);
        }
        c = parent.getComponent(3);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 312, 416, 24);
        }
        c = parent.getComponent(4);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 472, 96, 32);
        }
        c = parent.getComponent(5);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 432, 96, 32);
        }
        c = parent.getComponent(6);
        if (c.isVisible()) {
            c.setBounds(insets.left + 120, insets.top + 472, 328, 32);
        }
        c = parent.getComponent(7);
        if (c.isVisible()) {
            c.setBounds(insets.left + 424, insets.top + 152, 216, 24);
        }
        c = parent.getComponent(8);
        if (c.isVisible()) {
            c.setBounds(insets.left + 424, insets.top + 192, 216, 24);
        }
        c = parent.getComponent(9);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 232, 416, 24);
        }
        c = parent.getComponent(10);
        if (c.isVisible()) {
            c.setBounds(insets.left + 424, insets.top + 232, 216, 24);
        }
        c = parent.getComponent(11);
        if (c.isVisible()) {
            c.setBounds(insets.left + 424, insets.top + 312, 216, 24);
        }
        c = parent.getComponent(12);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 272, 416, 24);
        }
        c = parent.getComponent(13);
        if (c.isVisible()) {
            c.setBounds(insets.left + 424, insets.top + 272, 216, 24);
        }
        c = parent.getComponent(14);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 360, 144, 24);
        }
        c = parent.getComponent(15);
        if (c.isVisible()) {
            c.setBounds(insets.left + 160, insets.top + 360, 400, 24);
        }
        c = parent.getComponent(16);
        if (c.isVisible()) {
            c.setBounds(insets.left + 120, insets.top + 432, 136, 32);
        }
        c = parent.getComponent(17);
        if (c.isVisible()) {
            c.setBounds(insets.left + 544, insets.top + 472, 96, 32);
        }
        c = parent.getComponent(18);
        if (c.isVisible()) {
            c.setBounds(insets.left + 424, insets.top + 48, 216, 24);
        }
        c = parent.getComponent(19);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 48, 416, 24);
        }
        c = parent.getComponent(20);
        if (c.isVisible()) {
            c.setBounds(insets.left + 8, insets.top + 16, 616, 16);
        }
    }
}
