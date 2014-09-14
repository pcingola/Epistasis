package meshi.geometry;

import meshi.util.*;
import meshi.util.filters.Filter;
import java.util.*;

public class DistanceList implements Iterable<Distance>{
    protected  Distance[] internalArray;
    protected  int size;
    protected  int capacity;
    protected  Filter filter;


    public DistanceList(int capacity) {
	this(capacity,null);
    }
    public DistanceList(int capacity, Filter filter) {
	internalArray = new Distance[capacity];       
	size = 0;
	this.capacity = capacity;
	this.filter = filter;
    }
    
	
	
	
   public void clear() {
	if (capacity > size*1.01)
	    internalArray = new Distance[(int)1.01*size+1];
	capacity = internalArray.length;
	size = 0;
    }

    public boolean add(Distance distance) {
	if ((filter != null) && (! filter.accept(distance))) return false; 
 	if (size < internalArray.length) {
	    internalArray[size] = distance;
	    size++;
	    return true;
	}
	else {
	    capacity = 10+(capacity*10)/8;
	    Distance[] newArray = new Distance[capacity];
	    for (int i = 0; i < size; i++)
		newArray[i] = internalArray[i];
	    internalArray = newArray;
	    return add(distance);
	}
    }
    public Distance get(int index) {
	return internalArray[index];
    }
    
    public int size() {return size;}
    
    public boolean isEmpty() {return (size == 0);}
    public Distance[] toArray() {return internalArray;}

    public Iterator<Distance>iterator() {return new DistanceIterator(internalArray,size);}

    private static class DistanceIterator implements Iterator {
	private Distance[] array;
	private int current, size;
	public DistanceIterator(Distance[] array, int size) {
	    this.array = array;
	    this.size  = size;
	    current    = 0;
	}
	
	public Distance next() {
	    if (current >= size) throw new RuntimeException("no more distances");
	    current++;
	    return array[current-1];
	}

	public boolean hasNext() { return (current < size);}
	
	public void remove() {throw new RuntimeException("Do not do that");}
    }

    public void sort() {
	Arrays.sort(internalArray);
    }
	
}

