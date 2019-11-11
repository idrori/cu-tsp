/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.dataStructures;

public class DynamicArray implements Array {
	private Object[] arr;

	public DynamicArray(){
		arr = new Object[0];
	}

	public int size() {return Integer.MAX_VALUE;}

	public Object get(int i) {
		if (i < arr.length) return arr[i];
		else return null;
	}

public void set(int i, Object data) {
		ensureCapacity(i+1);
		arr[i] = data;
	}

	public void ensureCapacity(int capacity){
		if (capacity > arr.length){
			Object[] tmpArr = new Object[capacity];
			for (int j=0; j<arr.length; j=j+1){
				tmpArr[j] = arr[j];
			}
			arr = tmpArr;
		}
	}
} //
