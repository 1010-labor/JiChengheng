package mapreduce.wordcount;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.conf.Configured;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.IntWritable;
import org.apache.hadoop.io.RawComparator;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapred.JobConf;
import org.apache.hadoop.mapred.jobcontrol.*;
import org.apache.hadoop.mapred.jobcontrol.Job;
import org.apache.hadoop.mapreduce.*;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.input.TextInputFormat;
import org.apache.hadoop.mapreduce.lib.jobcontrol.JobControl;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.mapreduce.lib.output.TextOutputFormat;
import org.apache.hadoop.util.Tool;
import org.apache.hadoop.util.ToolRunner;

import java.net.URI;
import java.util.Map;

/**
 * Created by Chengheng on 2016/10/17.
 */
public class WordCountJob extends Configured implements Tool {
    public int run(String[] args) throws Exception {
        System.out.println("start");
        Configuration conf = getConf();
        org.apache.hadoop.mapreduce.Job job = org.apache.hadoop.mapreduce.Job.getInstance();
        System.out.println();
        Job cjob = new Job(new JobConf(conf));

//       job.setCombinerClass(WordConterReducer.class);
        //不是所有的mr都能用reducer代替Combiner： mean(a,b,c,d)!=mean(mean(a,b),mean(c,d))
        job.setMapperClass(WordCountMapper.class);
        job.setReducerClass(WordConterReducer.class);
        job.setNumReduceTasks(2);
        //两个reducer
        job.setMapOutputKeyClass(Text.class);
        job.setMapOutputValueClass(IntWritable.class);
        job.setInputFormatClass(TextInputFormat.class);  //设置输入方式
        job.setOutputFormatClass(TextOutputFormat.class);   //设置输出方式
        FileInputFormat.addInputPath(job, new Path(args[0]));
        FileOutputFormat.setOutputPath(job, new Path(args[1]));
        job.setJarByClass(WordCountJob.class);
        job.setPartitionerClass(WordCountPartitioner.class); //根据map输出key，分到不同的reducer上去，这里有两个reducer，所以WordCountPartitioner根据是否是字母开头的单词吧计数分到不同reducer
        job.setPriority(JobPriority.HIGH); //作业优先级
        job.setUser("datasafe");
        job.setSortComparatorClass(Text.Comparator.class); //map端排序，reduce端复制后的归并排序所用到的排序器
        cjob.setJob(job);
        int ret = job.waitForCompletion(true) ? 0 : 1;
        System.out.println("-------------------");
        System.out.println(ret);
        System.out.println(job.getCluster());
        System.out.println(job.getCounters());
        System.out.println(job.getFinishTime());
        System.out.println(job.getHistoryUrl());
        System.out.println(job.getJobFile());
        System.out.println(job.getJobName());
        System.out.println(job.getJobState());
        System.out.println(job.getPriority());
        System.out.println(job.getSchedulingInfo());
        System.out.println(job.getStartTime());
        System.out.println(job.getTrackingURL());
        System.out.println(job.getCacheFiles());
        System.out.println(job.getUser());
        System.out.println(job.getWorkingDirectory());
        return 0;
    }


}

