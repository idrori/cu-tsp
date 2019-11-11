/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.dataStructures;

public class FixedSizeArray implements Array {
	private Object[] arr;

	public FixedSizeArray(int size){
		arr = new Object[size];
	}
	public int size() {return arr.length ;}
	public Object get(int i) {return arr[i];}
	public void set(int i, Object data) {
		arr[i] = data;
	}
} //class FixedSizeArray
