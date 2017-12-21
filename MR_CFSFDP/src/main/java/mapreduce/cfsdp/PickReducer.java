package mapreduce.cfsdp;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer;

import java.io.IOException;
import java.util.*;

import static org.apache.hadoop.yarn.webapp.hamlet.HamletSpec.LinkType.index;

/**
 * Created by Chengheng on 2016/10/21.
 */

public class PickReducer extends Reducer<LongWritable, Text, LongWritable, Text> {
    static final int K_Density = 150;
    static final int K_NEIGHBOR = 20;
    static int ppp = 0;
    static int m = 0;

    @Override
    protected void reduce(LongWritable key, Iterable<Text> values, Context context) throws IOException, InterruptedException {
        int num = 0;
        List<Node> nodes = Lists.newArrayList();
        Map<Integer, Node> nodesMap = Maps.newHashMap();
        Map<Integer, Double> densityMap = Maps.newHashMap();
        Map<Integer, List<Integer>> neighborMap = Maps.newHashMap();
        for (Text value : values) {
            num++;
            String[] info = value.toString().split(" ");
            long gridNum = key.get();
            int index = Integer.valueOf(info[0]);
            double x = Double.valueOf(info[1]);
            double y = Double.valueOf(info[2]);
            boolean isBorder = new Boolean(info[3]);
            if (isBorder){
                System.out.println();
            }
            Node node = new Node(index, gridNum, x, y, isBorder);
            nodes.add(node);
            nodesMap.put(index, node);
        }
        for (int i = 0; i < nodes.size(); i++) {
            Node cur = nodes.get(i);
            if (cur.isborder) {
                m++;
                continue;
            }

            List<Node.Dist> dists = Lists.newArrayListWithCapacity(200);
            for (int j = 0; j < nodes.size(); j++) {
                if (i != j) dists.add(new Node.Dist(cur, nodes.get(j)));
            }
            double density = 0;
            Iterator<Node.Dist> it = dists.iterator();
            List<Integer> neighbors = Lists.newArrayListWithCapacity(K_NEIGHBOR);
            if (K_Density < num - 1) {
                //num-1个距离（不包含自身与自身的距离）
                Node.findKnn(dists, 0, num - 2, K_Density);
                Node.sortDist(dists, 0, K_Density);
                for (int k = 0; k < K_Density; k++) {
                    Node.Dist dist = it.next();
                    neighbors.add(dist.with);
                    density += dist.dist;
                }
                density = K_Density / density;
            } else {
                Node.sortDist(dists, 0, num - 2);
                for (int k = 0; k < num - 1; k++) {
                    Node.Dist dist = it.next();
                    if (k < K_NEIGHBOR) neighbors.add(dist.with);
                    density += dist.dist;
                }
                density = num == 1 ? 0 : num / density;
            }

            densityMap.put(cur.index, density);
            neighborMap.put(cur.index, neighbors);
        }
        for (Node node : nodes) {
            if (node.isborder) continue;
            Integer curIndex = node.index;
            double curDensity = densityMap.get(curIndex);
            boolean findHeavierNeiborInCurrentGrid = false;
            for (Integer neighborIndex : neighborMap.get(curIndex)) {
                if (densityMap.containsKey(neighborIndex) && densityMap.get(neighborIndex) > curDensity) {
                    findHeavierNeiborInCurrentGrid = true;
                    context.write(key, new Text("A " + curIndex + " " + densityMap.get(curIndex) + " " + neighborIndex));
                    break;
                }
            }
            if (!findHeavierNeiborInCurrentGrid) {
                Node cur = nodesMap.get(curIndex);
                context.write(key, new Text("B " + cur.index + " " + densityMap.get(curIndex) + " " + cur.x + " " + cur.y ));
            }


        }

    }


    @Override
    protected void cleanup(Context context) throws IOException, InterruptedException {
        System.out.println("-----------------m " + m);
        super.cleanup(context);
    }
}
