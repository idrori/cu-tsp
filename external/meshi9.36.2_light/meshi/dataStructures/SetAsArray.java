/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.dataStructures;

public class SetAsArray implements Set {
	private Array arr;
	int size;

	public SetAsArray(){
        this( new DynamicArray());
    }

    public SetAsArray(Array arr){
		this.arr = arr;
        if (arr.size() != 0)
            throw new RuntimeException("Array argument must be empty."+arr.size());
		size = 0;
	}

	public void add(Object data) {
		if (!contains(data)){
			arr.set(size, data);
			size = size+1;
		}
	}
    public void remove(Object data) {
            int i = indexOf(data);
            if (i>=0){
                arr.set(i, arr.get(size-1));
                arr.set(size-1, null);
                size = size-1;
            }
        }

        public boolean contains(Object data) {
            return indexOf(data) >= 0;
        }

    	public int size() {return size;}
    private int indexOf(Object data) {
            if (data != null){
                for (int i=0; i<size; i=i+1){
                    if (arr.get(i).equals(data)) return i;
                }
            }
            return -1;
        }
    } //class SetAsArray

