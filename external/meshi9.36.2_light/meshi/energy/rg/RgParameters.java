/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.rg;

/**
 *
 */
public class RgParameters {
    public final double alpha, beta, upperLimit, lowerLimit;
    public final String title;
    public static final RgParameters n_RGhSS = new RgParameters("log Number of atoms Vs. log RG of hydrophobic side chains in SS",
            4.356310e-001, -5.514293e-001, 1.097365e+000, 9.527014e-001);

//outliers 0 2 472 102
    public static final RgParameters RGhSS_RGhCoilCentered = new RgParameters("Log RG of hydrophobic side chains in SS Vs. log RG of hydrophobic side chains in Coil",
            5.805398e-001, 1.336618e+000, 1.141676e+000, 8.912878e-001);
//outliers 234 140 409 353
    public static final RgParameters RGhSS_RGbSScentered = new RgParameters("Log RG of hydrophobic side chains in SS Vs. log RG of CAs in SS",
            8.952122e-001, 3.613712e-001, 1.028366e+000, 9.718642e-001);
//outliers 799 477 250 41
    public static final RgParameters RGhSS_RGbCoilCentered = new RgParameters("Log RG of hydrophobic side chains in SS Vs. log RG of CAs in Coil",
            5.812256e-001, 1.395555e+000, 1.100033e+000, 9.141439e-001);
//outliers 364 384 155 257
    public static final RgParameters RGhSS_RGcSScentered = new RgParameters("Log RG of hydrophobic side chains in SS Vs. log RG of polar side chain getAtoms in SS",
            6.959428e-001, 1.087288e+000, 1.058440e+000, 9.263804e-001);
//outliers 104 237 278 309
    public static final RgParameters RGhSS_RGcCoilCentered = new RgParameters("Log RG of hydrophobic side chains in SS Vs. log RG of polar side chain getAtoms in Coil",
            5.342195e-001, 1.619181e+000, 1.110168e+000, 8.953153e-001);
//outliers 339 110 27 405
    public static final RgParameters RGhSS_RGbSS = new RgParameters("Log RG of hydrophobic side chains in SS Vs. log RG of backbone CAs getAtoms in SS", 9.029543e-001, 3.355479e-001, 1.028234e+000, 9.709244e-001);
//outliers 297 799 627 615
    public static final RgParameters RGhSS_RGcSS = new RgParameters("Log RG of hydrophobic side chains in SS Vs. log RG of polar side-chain getAtoms in SS", 7.321219e-001, 9.684661e-001, 1.057348e+000, 9.301641e-001);
//outliers 289 357 98 421

    public RgParameters(String title, double alpha, double beta, double upperLimit, double lowerLimit) {
        this.alpha = alpha;
        this.beta = beta;
        this.lowerLimit = lowerLimit;
        this.upperLimit = upperLimit;
        this.title = title;
    }
}
