/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.molecularElements.atoms.Atom;
import meshi.util.MeshiProgram;

public class Coordinates {
	public final static double PI = 3.14159265358979323846264338327950288;
	public static double NOWHERE_CONST = -9999.999;
	protected double[] x, y, z;
	private boolean random;
	private boolean nowhere;

	public Coordinates() {
		this(NOWHERE_CONST, NOWHERE_CONST, NOWHERE_CONST);
		nowhere = true;
	}

	public Coordinates(Atom atom) {
		this();
		if ((atom != null) && (!atom.nowhere())) {
			this.x[0] = atom.x();
			this.y[0] = atom.y();
			this.z[0] = atom.z();
			nowhere = atom.nowhere();
		}
	}

	public Coordinates(double x, double y, double z) {
		this.x = new double[2];
		this.y = new double[2];
		this.z = new double[2];
		this.x[0] = x;
		this.y[0] = y;
		this.z[0] = z;
		this.x[1] = 0.0;
		this.y[1] = 0.0;
		this.z[1] = 0.0;
		random = false;
        if (x == NOWHERE_CONST) {
            if ((y == NOWHERE_CONST) && (z == NOWHERE_CONST))
                nowhere = true;
            else throw new RuntimeException("This is weird: "+x+"  "+y+" "+z);
        }
        else
            nowhere = false;
	}

	public Coordinates(Coordinates coordinates) {
		this(coordinates.x(), coordinates.y(), coordinates.z());
		nowhere = coordinates.nowhere;
	}

	public Coordinates(Coordinates coordinates, double radius) {
		this(coordinates.x(), coordinates.y(), coordinates.z());
		double xx = radius, yy = radius, zz = radius;
		while (xx * xx + yy * yy + zz * zz > radius * radius) {
			xx = (2 * MeshiProgram.randomNumberGenerator().nextDouble() - 1)
					* radius;
			yy = (2 * MeshiProgram.randomNumberGenerator().nextDouble() - 1)
					* radius;
			zz = (2 * MeshiProgram.randomNumberGenerator().nextDouble() - 1)
					* radius;
		}
		addToX(xx);
		addToY(yy);
		addToZ(zz);
		random = true;
		nowhere = false;
	}

	public boolean random() {
		return random;
	}

	public final boolean nowhere() {
		return nowhere;
	}

	public void somewhere() {
		nowhere = false;
	}

	public void reset() {
		set(new Coordinates());
	}

	public final double x() {
		return x[0];
	}

	public final double y() {
		return y[0];
	}

	public final double z() {
		return z[0];
	}

	public final double fx() {
		return x[1];
	}

	public final double fy() {
		return y[1];
	}

	public final double fz() {
		return z[1];
	}

	// public void setX(double x) {this.x[0] = x;}

	// public void setY(double y) {this.y[0] = y;}
	// public void setZ(double z) {this.z[0] = z;}
	public void setXYZ(double x, double y, double z) {
		this.x[0] = x;
		this.y[0] = y;
		this.z[0] = z;
        nowhere = false;
	}

	public void set(Coordinates other) {
		x[0] = other.x[0];
		y[0] = other.y[0];
		z[0] = other.z[0];
		nowhere = other.nowhere;
	}

	public void set(Coordinates other, double radius) {
		set(new Coordinates(other, radius));
	}

	public void setFx(double fx) {
		this.x[1] = fx;
	}

	public void setFy(double fy) {
		this.y[1] = fy;
	}

	public void setFz(double fz) {
		this.z[1] = fz;
	}

	public void addToX(double addMe) {
		x[0] += addMe;
	}

	public void addToY(double addMe) {
		y[0] += addMe;
	}

	public void addToZ(double addMe) {
		z[0] += addMe;
	}

	public final void addToFx(double addMe) {
		x[1] += addMe;
	}

	public final void addToFy(double addMe) {
		y[1] += addMe;
	}

	public final void addToFz(double addMe) {
		z[1] += addMe;
	}

	public void resetForces() {
		x[1] = 0.0;
		y[1] = 0.0;
		z[1] = 0.0;
	}

	public final double[] X() {
		return x;
	}

	public final double[] Y() {
		return y;
	}

	public final double[] Z() {
		return z;
	}

	public boolean equals(Object obj) {
		Coordinates other = (Coordinates) obj;
		return ((this.z() == other.z()) & (this.y() == other.y()) & (this.x() == other
				.x()));
	}

	public String toString() {
		return "Coordinates (" + x() + "," + fx() + ")(" + y() + "," + fy()
				+ ")(" + z() + "," + fz() + ")";
	}

	public final double distanceFrom(Coordinates other) {
		double dx = x[0] - other.x[0];
		double dy = y[0] - other.y[0];
		double dz = z[0] - other.z[0];
		return Math.sqrt(dx * dx + dy * dy + dz * dz);
	}

	public final Coordinates MiddleCoordinates(Coordinates other) {
		double mx = (this.x() + other.x()) / 2;
		double my = (this.y() + other.y()) / 2;
		double mz = (this.z() + other.z()) / 2;
		Coordinates middCoordinates = new Coordinates(mx, my, mz);
		return middCoordinates;
	}

	public final Coordinates differenceVector(Coordinates other) {
		double dx = other.x() - this.x();
		double dy = other.y() - this.y();
		double dz = other.z() - this.z();
		Coordinates diffCoordinates = new Coordinates(dx, dy, dz);
		return diffCoordinates;
	}

	// Additional methods for vector calculus, Tommer 30.11.14

	public double norm() {
		return Math.sqrt(this.x() * this.x() + this.y() * this.y() + this.z()
				* this.z());
	}

	public final Coordinates vectorAdd(Coordinates vec) {
		return new Coordinates(this.x() + vec.x(), this.y() + vec.y(), this.z()
				+ vec.z());
	}

	public final Coordinates scalarMult(double a) {
		return new Coordinates(this.x() * a, this.y() * a, this.z() * a);
	}

	public final Coordinates scalarDivision(double a) {
		if (a == 0)
			throw new RuntimeException("can't divide by zero");
		return new Coordinates(this.x() / a, this.y() / a, this.z() / a);
	}

	public final Coordinates normalize() {
		return new Coordinates(this.x() / this.norm(), this.y() / this.norm(),
				this.z() / this.norm());
	}

	public final Coordinates vectorMult(Coordinates vec) {
		return new Coordinates(this.y() * vec.z() - this.z() * vec.y(),
				this.z() * vec.x() - this.x() * vec.z(), this.x() * vec.y()
						- this.y() * vec.x());
	}

	public static double sineOfOut_of_PlaneAngle(double angle1st_to_target,
			double angle2nd_to_target, double middleAngle) {
		if (Math.sin(middleAngle) == 0)
			throw new RuntimeException(
					"sin(middleAngle) is zero. this shouldn't happen");
		if ((2 * PI - angle1st_to_target - angle2nd_to_target <= middleAngle)||
				(-angle1st_to_target+middleAngle<=-angle2nd_to_target)||
				(-angle2nd_to_target+middleAngle<=-angle1st_to_target)) {
			return 0;// this is to take care of cases were the getAtoms are in one
						// plane without risk
			// of having to take the root of a negative quantity.
			// it is assumed that the angles are given so that middleAngle
			// cannot be significantly
			// larger than the quantity inside the if statement. this is a
			// geometric necessity.
		}
		return 1/ Math.sin(middleAngle)
				* Math.sqrt(1 - Math.cos(angle1st_to_target)
						* Math.cos(angle1st_to_target)
						- Math.cos(angle2nd_to_target)
						* Math.cos(angle2nd_to_target) - Math.cos(middleAngle)
						* Math.cos(middleAngle) + 2
						* Math.cos(angle1st_to_target)
						* Math.cos(angle2nd_to_target) * Math.cos(middleAngle));
	}

	public static double cosineOfgamma2(//needs more explaining! Tommer 24.1.15
			double sineTheta, double angle1st_to_target) {
		return -Math.cos(angle1st_to_target)
				/ Math.sqrt(1 - sineTheta * sineTheta);
		// here sineTheta should be given by the result of the previous function:
		//sineOfOut_of_PlaneAngle
		// and angle1st_to_target should be the angle between one of the lines
		// connecting an atom with
		// known coordinates to a middle atom (preferably also with known
		// coordinates) and the line
		// connecting the middle atom to an atom with unknown coordinates, which
		// we are trying to locate.
		// note: we have that angle from angleParameters.dat
	}

	

}


