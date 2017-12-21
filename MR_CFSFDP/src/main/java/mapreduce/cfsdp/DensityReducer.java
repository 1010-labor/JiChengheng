/*
package mapreduce.cfsdp;

import com.google.common.collect.Lists;
import com.google.common.util.concurrent.AtomicDouble;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Counter;
import org.apache.hadoop.mapreduce.ReduceContext;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.task.reduce.ExceptionReporter;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

*/
/**
 * Created by Chengheng on 2016/10/21.
 *//*

public class DensityReducer extends Reducer<IntWritable, Text, IntWritable, Text> {
    static final int K_Density = 30;
    static final int K_NEIGHBOR = 20;
    double max_density = 0;
    static AtomicDouble global_max_density = new AtomicDouble(0);

    @Override
    protected void reduce(IntWritable key, Iterable<Text> values, Context context) throws IOException, InterruptedException {
        List<Node> nodes = Lists.newArrayList();
        int num = -1;
        Counter NUM_OF_THIS_Reucer = context.getCounter("NUM_OF_POINTS", key.toString());
        for (Text value : values) {
            num++;

            try {
                String[] info = value.toString().split(" ");
                String isBorder = info[3];
                int index = Integer.valueOf(info[0]);
                int gridNum = key.get();
                double x = Double.valueOf(info[1]);
                double y = Double.valueOf(info[2]);

                Node node = new Node(index, gridNum, x, y, isBorder);
                nodes.add(node);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(-1);
            }

        }
        NUM_OF_THIS_Reucer.setValue(num + 1);
        Iterator<Node> it1 = nodes.iterator();
        while (it1.hasNext()) {
            Node cur = it1.next();
            if (cur.isborder) continue;
            List<Node.Dist> dists = Lists.newArrayListWithCapacity(1000);
            Iterator<Node> it2 = nodes.iterator();
            while (it2.hasNext()) {
                dists.add(new Node.Dist(cur, it2.next()));
            }
            double density = 0;
            Iterator<Node.Dist> it = dists.iterator();
            StringBuilder append = new StringBuilder();
            append.append(cur.gridnum);
            append.append(" " + cur.isborder);
            append.append(" " + cur.index);
            append.append(" " + cur.x);
            append.append(" " + cur.y);


            if (K_Density < num) {
                Node.findKnn(dists, 0, num, K_Density);
                Node.sortDist(dists, 0, K_Density);
                for (int i = 0; i < K_Density; i++) {
                    Node.Dist dist = it.next();
                    density += dist.dist;
                    if (i < K_NEIGHBOR) append.append(" " + dist.with);
                }
                append.append(" " + density / K_Density);
            } else {
                Node.sortDist(dists, 0, num);
                for (int i = 0; i <= num; i++) {
                    Node.Dist dist = it.next();
                    density += dist.dist;
                    if (i < K_NEIGHBOR) append.append(" " + dist.with);
                }
                append.append(" " + density / (num + 1));
            }
            context.write(key,new Text(append.toString()));
            if (this.max_density < cur.density) {
                this.max_density = cur.density;
            }
        }

    }



    @Override
    protected void cleanup(Context context) throws IOException, InterruptedException {
        super.cleanup(context);
    }
}
*/
