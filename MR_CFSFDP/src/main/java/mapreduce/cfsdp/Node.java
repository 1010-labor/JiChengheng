package mapreduce.cfsdp;

import org.apache.commons.math3.analysis.function.Pow;
import org.apache.commons.math3.analysis.function.Sqrt;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by Chengheng on 2016/10/21.
 */
public class Node {
    public int index;
    public long gridnum;
    public double x;
    public double y;
    boolean isborder;
    public double density;
    public List<Integer> neighbors;
    public double deta;
    public double gama;

    Node(int index, long gridnum, double x, double y, boolean isborder) {
        this.index = index;
        this.gridnum = gridnum;
        this.x = x;
        this.y = y;
        this.isborder = isborder;
    }

    public static class Dist implements Comparable<Dist> {
        int with;
        double dist;
        static Pow pow = new Pow();
        static Sqrt sqrt = new Sqrt();

        Dist(Node a, Node b) {
            this.with = b.index;
            this.dist = sqrt.value(pow.value((a.x - b.x), 2) + pow.value((a.y - b.y), 2));
        }

        Dist(int with, double dist) {
            this.dist = dist;
            this.with = with;
        }

        public int compareTo(Dist o) {
            return this.dist > o.dist ? 1 : 0;
        }
    }

    static public void findKnn(List<Dist> dm_line, int left, int right, int k) {

        int current_find = right - left + 1;
        Dist key;
        int lindex;
        int rindex = right;
        int key_history = 0;

        while (current_find > k + 1) {

            lindex = left;
            key_history = rindex;
            key = dm_line.get(left);
            while (lindex < rindex) {
                while (lindex < rindex && key.dist <= dm_line.get(rindex).dist) rindex--;
                if (lindex < rindex) dm_line.set(lindex++, dm_line.get(rindex));
                while (lindex < rindex && key.dist >= dm_line.get(lindex).dist) lindex++;
                if (lindex < rindex) dm_line.set(rindex--, dm_line.get(lindex));
            }
            assert (lindex == rindex);
            dm_line.set(lindex, key);
            current_find = rindex - left + 1;
        }
        if (current_find == k + 1 || current_find == k) return;
        else {
            findKnn(dm_line, rindex + 1, key_history, k - current_find);
        }
    }

    static public void sortDist(List<Dist> dm_line, int left, int right) {
        if (right > left) {
            Dist key = dm_line.get(left);
            int lindex = left;
            int rindex = right;
            while (lindex < rindex) {
                while (lindex < rindex && key.dist <= dm_line.get(rindex).dist) rindex--;
                if (lindex < rindex) dm_line.set(lindex++, dm_line.get(rindex));
                while (lindex < rindex && key.dist >= dm_line.get(lindex).dist) lindex++;
                if (lindex < rindex) dm_line.set(rindex--, dm_line.get(lindex));
            }
            dm_line.set(lindex, key);
            sortDist(dm_line, left, lindex - 1);
            sortDist(dm_line, rindex + 1, right);
        } else {
            return;
        }

    }

    public static void main(String[] args) {
        List<Dist> l = new ArrayList<Dist>();
    }
}
