ACM-cheat-sheet
===============

#基础 Basic
####Buffered Input
```java
import java.io.*;
```
```java
BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
    int T = Integer.parseInt(br.readLine());
    String[] data;
    while (T-- > 0) {
        data = br.readLine().split(" ");
        n = Integer.parseInt(data[0]);
        m = Integer.parseInt(data[1]);
```

####Integer & BigInt
```java
```

####建图
```java
class Edge {
    public int e;
    public int l;
    public Edge(int e, int l) {
        this.e = e;
        this.l = l;
    }
}
public static ArrayList<ArrayList<Edge>> map;
```
```java
for (int i = 0; i < N; i++) {
    ArrayList<Edge> newStart = new ArrayList<Edge>();
    map.add(newStart);
}
```
```java
for (int i = 0; i < M; i++) {
    data = br.readLine().split(" ");
    s = Integer.parseInt(data[0]);
    e = Integer.parseInt(data[1]);
    l = Integer.parseInt(data[2]);
    Edge edge = new Edge(e - 1, l); map.get(s-1).add(edge);
    edge = new Edge(s - 1, l); map.get(e-1).add(edge);
}
```

#数据结构 Data Structure

####TreeSet
```java
//Constructor and Description
TreeSet()
Constructs a new, empty tree set, sorted according to the natural ordering of its elements.

//Modifier and Type	Method and Description
boolean	add(E e)
Adds the specified element to this set if it is not already present.
boolean	contains(Object o)
Returns true if this set contains the specified element.
Iterator<E>	descendingIterator()
Returns an iterator over the elements in this set in descending order.
boolean	isEmpty()
Returns true if this set contains no elements.
Iterator<E>	iterator()
Returns an iterator over the elements in this set in ascending order.
boolean	remove(Object o)
Removes the specified element from this set if it is present.
int	size()
Returns the number of elements in this set (its cardinality).
```


####TreeMap
```java
//Constructor and Description
TreeMap()
Constructs a new, empty tree map, using the natural ordering of its keys.

//Modifier and Type    Method and Description
boolean	containsKey(Object key)
Returns true if this map contains a mapping for the specified key.

boolean	containsValue(Object value)
Returns true if this map maps one or more keys to the specified value.

K	firstKey()
Returns the first (lowest) key currently in this map.

V	get(Object key)
Returns the value to which the specified key is mapped, or null if this map contains no mapping for the key.

Set<K>	keySet()
Returns a Set view of the keys contained in this map.

K	lastKey()
Returns the last (highest) key currently in this map.

V	put(K key, V value)
Associates the specified value with the specified key in this map.

void	putAll(Map<? extends K,? extends V> map)
Copies all of the mappings from the specified map to this map.

V	remove(Object key)
Removes the mapping for this key from this TreeMap if present.

int	size()
Returns the number of key-value mappings in this map.
```

####线段树(Segment Tree)
```java
void initialize(intnode, int b, int e, int M[MAXIND], int A[MAXN], int N){
    if (b == e)
        M[node] = b;
    else {
        //compute the values in the left and right subtrees
        initialize(2 * node, b, (b + e) / 2, M, A, N);
        initialize(2 * node + 1, (b + e) / 2 + 1, e, M, A, N);
        //search for the minimum value in the first and 
        //second half of the interval
        if (A[M[2 * node]] <= A[M[2 * node + 1]])
            M[node] = M[2 * node];
        else
            M[node] = M[2 * node + 1]; 
    }
} 
```
```java
int query(int node, int b, int e, int M[MAXIND], int A[MAXN], int i, int j) {
    int p1, p2;
    //if the current interval doesn't intersect 
    //the query interval return -1
    if (i > e || j < b)
        return -1;
    //if the current interval is included in 
    //the query interval return M[node]
    if (b >= i && e <= j)
        return M[node];
   
    //compute the minimum position in the 
    //left and right part of the interval
    p1 = query(2 * node, b, (b + e) / 2, M, A, i, j);
    p2 = query(2 * node + 1, (b + e) / 2 + 1, e, M, A, i, j);
    //return the position where the overall 
    //minimum is
    if (p1 == -1)
        return M[node] = p2;
    if (p2 == -1)
        return M[node] = p1;
    if (A[p1] <= A[p2])
        return M[node] = p1;
    return M[node] = p2;
}
```
####RMQ(Range Minimum Query)
[for LCA(Lowest Common Ancestor)](http://community.topcoder.com/tc?module=Static&d1=tutorials&d2=lowestCommonAncestor)
```java
void process2(int M[MAXN][LOGMAXN], int A[MAXN], int N) {
    int i, j;
    //initialize M for the intervals with length 1
    for (i = 0; i < N; i++)
        M[i][0] = i;
    //compute values from smaller to bigger intervals
    for (j = 1; 1 << j <= N; j++)
        for (i = 0; i + (1 << j) - 1 < N; i++)
            if (A[M[i][j - 1]] < A[M[i + (1 << (j - 1))][j - 1]])
                M[i][j] = M[i][j - 1];
            else
                M[i][j] = M[i + (1 << (j - 1))][j - 1];
}


public static int rmq(int l, int r) {
    int j = (int)(Math.log(r - l + 1) / Math.log(2));
    int len = 1 << j;
    return nodeDepth[dp[l][j]] < nodeDepth[dp[r - len + 1][j]] ? dp[l][j] : dp[r - len + 1][j];
}

public static int lca(int x, int y) {
    int ans;
    if (nodePosition[x] > nodePosition[y]) {
        ans = rmq(nodePosition[y], nodePosition[x]);
    } else {
        ans = rmq(nodePosition[x], nodePosition[y]);
    }
    return visitOder[ans] + 1;
}
```
```java
public static void dfs(int currentNode, int depth) {
        visited[currentNode] = true;
        nodeDepth[tot] = depth;
        visitOder[tot] = currentNode;
        nodePosition[currentNode] = tot++;
        for (int child : map.get(currentNode)) {
            if (!visited[child]) {
                dfs(child, depth + 1);
                visitOder[tot] = currentNode;
                nodeDepth[tot++] = depth;
            }
        }
    }
```

#算法 algorithm
####快排 (quickSort)
```java
public static void qsort(int l, int r) {
    int i = l;
    int j = r;
    int x = edge[(l + r) / 2];
    int tmp;
    while (true) {
        while (edge[i] > x) i++;
        while (edge[j] < x) j--;
        if (i <= j) {
            swap(i, j);
            i++;
            j--;
        }
        if (i > j)
            break;
    }
    if (j > l) qsort(l, j);
    if (i < r) qsort(i, r);
}
```
####差分约束系统(?)
        不会，需要学

####Astar搜索(A*搜索)
        还不会，需要学

####广搜 + 最短路 (SPFA)
```java
dist = new int[N + 1];
    f = new boolean[N + 1];
    time = new int[N + 1];
    for (int i = 0; i < N; i++) {
        dist[i] = 1147483648;
        f[i] = false;
        time[i] = 0;
    }
    Queue<Integer> q = new LinkedList<Integer>();
    q.add(N);
    f[N] = true;
    while (q.peek() != null) {
        start = q.poll();
        for (Edge end : map.get(start)){
            if (dist[start] + end.l < dist[end.e]) {
                dist[end.e] = dist[start] + end.l;
                if (!f[end.e]) {
                    q.add(end.e);
                    f[end.e] = true;
                }
                time[end.e]++;
                if (time[end.e] > N) {
                    System.out.println("YES");
                    return;
                }
            }
        }
        f[start] = false;
    }
```

####最大流 Maxflow
```
I will change the code into java asap
```
```pascal
program asdf;
  var
    x,y,n,m,z,i,flow,ans:longint;
    d,num,t:array[0..200] of longint;
    g,p:array[0..200,0..30] of longint;
  function min(a,b:longint):longint;
    begin
      if a>b then
        exit(b)
      else exit(a);
    end;
  function dfs(u,flow:longint):longint;
    var
      v,now,i:longint;
    begin
      if u=m then
        exit(flow);
      dfs:=0;
      for i:=1 to t[u] do
        begin
          v:=p[u,i];
          if (g[u,v]>0)and(d[u]=d[v]+1) then
            begin
              now:=dfs(v,min(flow-dfs,g[u,v]));
              dec(g[u,v],now);
              inc(g[v,u],now);
              dfs:=dfs+now;
              if dfs=flow then
                exit(dfs);
            end;
          end;
      if d[1]>=m then
        exit;
      dec(num[d[u]]);
      if num[d[u]]=0 then
        d[1]:=m;
      inc(d[u]);
      inc(num[d[u]]);
    end;
  begin
    while not eof do
      begin
        readln(n,m);
        fillchar(t,sizeof(t),0);
        fillchar(g,sizeof(g),0);
        for i:=1 to n do
          begin
            readln(x,y,z);
            g[x,y]:=g[x,y]+z;
            inc(t[x]);
            p[x,t[x]]:=y;
            inc(t[y]);
            p[y,t[y]]:=x;
          end;
        fillchar(d,sizeof(d),0);
        num[0]:=m;
        ans:=0;
        while d[1]<m do
          begin
            flow:=dfs(1,maxlongint);
            ans:=ans+flow;
          end;
        writeln(ans);
      end;
  end.
```

####最小费用最大流
```
最小费用最大流的算法
基本思路：
    把弧<i,j>的单位费用w[i,j]看作弧<i,j>的路径长度，每次找从源点s到汇点t长度最短（费用最小）的可增广路径进行增广。
1. 最小费用可增广路
2. 路径s到t的长度即单位流量的费用。
ps：是网络流EK算法的改进，在求增广路径的时候，把bfs改为带权的spfa，每次求权值最小的增广路。
ps：要注意一点，逆边cost[i][j] = -cost[j][i]，不能忘了加上去
```


#数学 Mathematic

####Miller-Rabbin Prime test
```java
public static boolean isPrime(int n) {
    int a;
    for (int i = 0; i < times; i++) {
        a = 2 + (int)(Math.random() * (n - 2));
        if (modular_exp(a, n-1, n) != 1) {
            return false;
        }
    }
    return true;
}
```
####快速幂 Fast Exponention
```java
public static int modular_exp(int a, int d, int n) {
    int ans = 1;
    while (d > 0) {
        if (d % 2 == 1) {
            ans = ans * a % n;
        }
        a = a * a % n;
        d /= 2;
    }
    return ans % n;
}
```
