package com.insilico.application.insilicotools.util;

import com.google.common.base.Function;

import java.util.Vector;
import java.util.concurrent.*;

/**
 * Created by jfeng1 on 10/17/15.
 */
public class FutureWorker {
    Vector<Future> futures;
    ExecutorService threadPool;

    public FutureWorker(int numThreads, final Vector<String> input, final Function function) {
        futures = new Vector<Future>();
        threadPool = Executors.newFixedThreadPool(numThreads);
        for(int i=0;i<numThreads;i++){
            final int idx = i;
            futures.add(threadPool.submit(new Callable() {
                @Override
                public Object call() throws Exception {
                    System.out.println(""+idx);
                    int sum = 0;
                    for(int i=0;i<100000/(idx+1);i++){
                        sum += i;
                    }
                    System.out.println(idx + " "+sum);
                    return ""+idx+" finished.";
                }
            }));
        }

    }

    public void runAll() throws ExecutionException, InterruptedException {
        for(Future task:futures){
            System.out.println("Running ...");
            if(task.isDone()){
                System.out.println((String)task.get());
            }else {
                Thread.sleep(1);
                System.out.println("idling ...");
            }
        }
        threadPool.shutdown();
    }



    public static void main(String[] args) {
        Vector<String> input = new Vector<String>();
        input.add("Future 1");
        input.add("Future 2");
        input.add("Future 3");
        input.add("Future 4");
        input.add("Future 5");

        FutureWorker w = new FutureWorker(5, input, new Function() {
            @Override
            public Object apply(Object o) {
                String i1 = (String)o;
                return i1+" is done.";
            }
        });
        try {
            w.runAll();

        } catch (ExecutionException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }
}
