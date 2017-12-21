package mapreduce.cfsdp;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import org.apache.commons.math3.analysis.function.Pow;
import org.apache.commons.math3.analysis.function.Sqrt;

import java.io.*;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

/**
 * Created by Chengheng on 2016/10/24.
 * 第二趟mapreduce的替代，从第一趟mapreduce的结果中取出代表样本，进行聚类，并为其他样本归类
 */
public class RES {
    public static void main(String[] args) throws IOException {
        BufferedReader bfr1 = new BufferedReader(new FileReader(new File("out/part-r-00000")));
        String line = bfr1.readLine();

        Map<Integer, Integer> nei = Maps.newHashMap();

        List<Node> all = Lists.newArrayList();
        List<Node> B_nodes = Lists.newArrayList();

        while (line != null) {
            String[] tokens = line.split(" ");
            Node nn = new Node(Integer.valueOf(tokens[1]), -1, 0, 0, false);
            nn.density = Double.valueOf(tokens[2]);
            all.add(nn);
            if (tokens.length > 4) {
                Node n = new Node(Integer.valueOf(tokens[1]), -1, Double.valueOf(tokens[3]), Double.valueOf(tokens[4]), false);
                n.density = Double.valueOf(tokens[2]);
                B_nodes.add(n);
            } else {
                nei.put(Integer.valueOf(tokens[1]), Integer.valueOf(tokens[3]));
            }
            line = bfr1.readLine();
        }
        Collections.sort(B_nodes, new Comparator<Node>() {
            public int compare(Node o1, Node o2) {
                return -1 * Double.compare(o1.density, o2.density);
            }
        });

        Pow pow = new Pow();
        Sqrt sqrt = new Sqrt();
        double max_deta = 0;
        double max_density = 0;
        double min_deta = 1000000;
        double min_density = 1000000;
        for (int i = 0; i < B_nodes.size(); i++) {
            Node thisNode = B_nodes.get(i);
            max_density = thisNode.density >= max_density ? thisNode.density : max_density;
            min_density = min_density <= thisNode.density ? min_density : thisNode.density;
            double deta = 1000000;
            for (int j = 0; j < i; j++) {
                Node neiborNode = B_nodes.get(j);
                double tmp = sqrt.value(pow.value((neiborNode.x - thisNode.x), 2) + pow.value((neiborNode.y - thisNode.y), 2));
                if (tmp < deta && neiborNode.density > thisNode.density) {
                    deta = tmp;
                    nei.put(thisNode.index, neiborNode.index);
                }
            }
            if (i > 0) {
                thisNode.deta = deta;
                max_deta = deta > max_deta ? deta : max_deta;
                min_deta = deta <= min_deta ? deta : min_deta;
            }
        }
        B_nodes.get(0).deta = max_deta;
        for (Node b_node : B_nodes) {
            b_node.gama = b_node.density * b_node.deta;
        }
        Collections.sort(B_nodes, new Comparator<Node>() {
            public int compare(Node o1, Node o2) {
                return -1 * Double.compare(o1.gama, o2.gama);
            }
        });



        int[] clus = new int[19200];


        for (int i = 0; i < 4; i++) {
            int t = B_nodes.get(i).index;
            clus[t] = i + 1;
        }
        Collections.sort(B_nodes, new Comparator<Node>() {
            public int compare(Node o1, Node o2) {
                return -1 * Double.compare(o1.density, o2.density);
            }
        });
        System.out.println("---------------");
        System.out.println(B_nodes.size());
        System.out.println("---------------");
        for (int i = 1; i < B_nodes.size(); i++) {
            if (clus[B_nodes.get(i).index] == 0) {
                int r = clus[nei.get(B_nodes.get(i).index)];
                clus[B_nodes.get(i).index] = r;
            }
        }

        Collections.sort(all, new Comparator<Node>() {
            public int compare(Node o1, Node o2) {
                return -1 * Double.compare(o1.density, o2.density);
            }
        });
       for (Node node : all) {
            if (clus[node.index] == 0) {




                 clus[node.index] = clus[nei.get(node.index)];
            }
        }
        FileWriter fw = new FileWriter(new File("Working/res"));
        for (int i = 0; i < 19200; i++) {
            fw.write(i + " " + clus[i] + " 0\n");
        }
        fw.close();
    }
}
