ACM-cheat-sheet
===============

#基础 Basic
####Buffered Input/Output
```java
import java.io.*;
```
```java
BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
BufferedWriter bw = new BufferedWriter(new OutputStreamReader(System.out));
	
int T = Integer.parseInt(br.readLine());
String[] data;
while (T-- > 0) {
	data = br.readLine().split(" ");
       	n = Integer.parseInt(data[0]);
       	m = Integer.parseInt(data[1]);
       	...
}

br.close();
bw.close();
```

####Integer & BigInt (very useful)
```java
java.math.BigInteger
//Constructor and Description
BigInteger(String val)
Translates the decimal String representation of a BigInteger into a BigInteger.
BigInteger(String val, int radix)
Translates the String representation of a BigInteger in the specified radix into a BigInteger.

//Modifier and Type	  Method and Description
BigInteger	abs()
Returns a BigInteger whose value is the absolute value of this BigInteger.

BigInteger	add(BigInteger val)
Returns a BigInteger whose value is (this + val).
BigInteger	subtract(BigInteger val)
BigInteger	multiply(BigInteger val)
BigInteger	divide(BigInteger val)
BigInteger	pow(int exponent)
Returns a BigInteger whose value is (thisexponent).

BigInteger	and(BigInteger val)
Returns a BigInteger whose value is (this & val).
BigInteger	or(BigInteger val)
BigInteger	not()

BigInteger	shiftLeft(int n)
Returns a BigInteger whose value is (this << n).
BigInteger	shiftRight(int n)

BigInteger	mod(BigInteger m)
Returns a BigInteger whose value is (this mod m).
BigInteger	remainder(BigInteger val)
Returns a BigInteger whose value is (this % val).
BigInteger	gcd(BigInteger val)
Returns a BigInteger whose value is the greatest common divisor of abs(this) and abs(val).

BigInteger	max(BigInteger val)
Returns the maximum of this BigInteger and val.
BigInteger	min(BigInteger val)

int	compareTo(BigInteger val)
Compares this BigInteger with the specified BigInteger.
String	toString(int radix)
Returns the String representation of this BigInteger in the given radix.
static BigInteger	valueOf(long val)
Returns a BigInteger whose value is equal to that of the specified long.
long	longValue()
Converts this BigInteger to a long.
```

####建图
```java
class Edge {
    public int e; public int l;
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
java.util.TreeSet<E>

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
java.util.TreeMap<K,V>
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
Returns the value to which the specified key is mapped
    or null if this map contains no mapping for the key.
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
####堆 (Heap)
```java
public static void buildheap() {
    n = arr.length-1;
    for(int i = n/2; i >= 0; i--){
        maxheap(i);
    }
}
public static void maxheap(int i) {
    left = 2 * i; right = 2 * i + 1;
    largest = i;
    if(left <= n && arr[left] > arr[i]) {
        largest = left;
    } else {
        if(right <= n && arr[right] > arr[largest])
            largest=right;
    }
    if(largest != i) {
        exchange(i, largest); //i and l are index
        maxheap(largest);
    }
}
public static void sort() {
    buildheap();
    for(int i = n; i > 0;i--) {
        exchange(0, i);
        n = n - 1;
        maxheap(0);
    }
}

```
####线段树(Segment Tree)
```java
void initialize(intnode, int b, int e, int M[MAXIND], int A[MAXN], int N) {
    if (b == e)
        M[node] = b;
    else {
        initialize(2 * node, b, (b + e) / 2, M, A, N);
        initialize(2 * node + 1, (b + e) / 2 + 1, e, M, A, N);
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
    //change the edge to your own array
}
```

####最小生成树 prime
```java
public int solve() {
	int i,j,Min,v,sum=0;
	//每次加入一个节点
	for(i=1;i<n;i++) {
		Min=MAX.num;
		v=0;
		for(j=1;j<=n;j++) 
		if(visit[j]==0&&dis[j]<Min) {
			Min=dis[j];
			v=j;
		}
		sum+=Min;
		visit[v]=1;
		for(j=1;j<=n;j++)
     		if(visit[j]==0&&dis[j]>map[v][j])
	   	dis[j]=map[v][j];
	}
	return sum;
}
```
####最小生成树 kruskal

####差分约束系统(?)
如果一个系统由n个变量和m个约束条件组成，其中每个约束条件形如xj-xi<=bk(i,j∈[1,n],k∈[1,m]),则其为差分约束系统(system of difference constraints)。亦即，差分约束系统是关于一组变量的特殊不等式组。求解差分约束系统，可以转化成图论的单源最短路径问题。
观察xj-xi<=bk，会发现它类似最短路中的三角不等式d[v]<=d[u]+w[u,v]，即d[v]-d[u]<=w[u,v]。因此，以每个变量xi为结点，对于约束条件xj-xi<=bk，连接一条边(i,j)，边权为bk。我们再增加一个源点s,s与所有点相连，边权均为0。对这个图，以s为源点运行bellman-ford算法，最终{d[i]}即为一组可行解。(差分约束系统的解的一个特点是，当将所有变量同时增加相同的大小，约束条件依然成立）

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
####Dijkstra最短路 Java可快速加入二叉堆优化？
```python
for i=1:n
     ins=0;
     for j=1:length(s)
        if i==s(j)
           ins=1;
     end,  end
     if ins==0
        v=i;
        if label(v)>(label(u)+w(u,v))
           label(v)=(label(u)+w(u,v)); f(v)=u;
  end, end, end   
```

####二分图 匈牙利算法
```c
int path(int u) {
	int v;
	for(v=0;v<ny;v++) {
		if(edge[u][v]&&!visited[v]) {
			visited[v]=1;
            if(cy[v]==-1||path(cy[v])) {
            //如果y集合中的v元素没有匹配或者是v已经匹配，但是从cy[v]中能够找到一条增广路
            	cx[u]=v;
            	cy[v]=u;
            	return 1;
            }
      	}
	}
	return 0;
}
=======
//main|
//=====
for(int i=0;i<=nx;i++) {
	if(cx[i]==-1) {
		memset(visited,0,sizeof(visitited));
		res+=path(i);
	}
}
```
####最大流 Maxflow
```
I will change the code into java asap
```
```pascal
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
            dec(g[u,v],now); inc(g[v,u],now);
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
          inc(t[x]); p[x,t[x]]:=y;
          inc(t[y]); p[y,t[y]]:=x;
        end;
      fillchar(d,sizeof(d),0);
      num[0]:=m; ans:=0;
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
```cpp
#define maxn 1005
#define inf 0x3f3f3f3f
struct edge{int v,w,f,c,next;} e[50000];
int vit[maxn],dis[maxn],start[maxn],p[maxn];
int tot,n,m;
void _add(int v,int w,int f,int c) {
    e[tot].v=v; e[tot].w=w; e[tot].f=f; e[tot].c=c;
    e[tot].next=start[v];start[v]=tot++;
}
void add(int v,int w,int f,int c) {
    _add(v,w,f,c);
    _add(w,v,0,-c);
}
bool spfa(int s,int t,int n) {
    int v,w;	//寻找费用最小的可增广路
    queue<int> q;
    for(int i=0;i<n;i++)
    { p[i]=-1;vit[i]=0;dis[i]=inf; }
    vit[s]=1;dis[s]=0;q.push(s);
    while(!q.empty()) {
        v=q.front();q.pop();vit[v]=0;
        for(int i=start[v];i!=-1;i=e[i].next)
            if(e[i].f) {
                w=e[i].w;
                if(dis[w]>dis[v]+e[i].c) {
                    dis[w]=dis[v]+e[i].c;
                    p[w]=i;
                    if(!vit[w]) {vit[w]=1;q.push(w);}
                }
            }
    }
    return dis[t]!=inf;
}
int cost(int s,int t,int n) {
    int ans=0,flow=inf,i;
    while(spfa(s,t,n)) {
        ans+=dis[t];
        for(i=p[t];i!=-1;i=p[e[i].v])//可改进量
            if(e[i].f<flow) flow=e[i].f;
        for(i=p[t];i!=-1;i=p[e[i].v]) {
            e[i].f-=flow;//调整
            e[i^1].f+=flow;
        }
    }
    return ans;
}
int main() {
    int i,v,w,c;
    scanf("%d%d",&n,&m);
    for(i=0;i<n+2;i++) start[i]=-1;tot=0;//初始化
    for(i=0;i<m;i++) {
        scanf("%d%d%d",&v,&w,&c);
        add(v,w,1,c);//添加边 此题为双向边
        add(w,v,1,c);
    }
    add(0,1,2,0);add(n,n+1,2,0);//此题添加的虚拟源点汇点
    printf("%d\n",cost(0,n+1,n+2));//最小费用
    //n+1 is the end of flow. n+2 is the number of nodes.
    return 0;
}
```

####最小割最大流

    最重要的在于建图
    一般题目会比较难

####强连通分量 Tarjan
```cpp
void tarjan(int i) {
    int j;
    DFN[i]=LOW[i]=++Dindex;
    instack[i]=true;
    Stap[++Stop]=i;
    for (edge *e=V[i];e;e=e->next) {
        j=e->t;
        if (!DFN[j]) {
            tarjan(j);
            if (LOW[j]<LOW[i])
                LOW[i]=LOW[j];
        }
        else if (instack[j] && DFN[j]<LOW[i])
            LOW[i]=DFN[j];
    }
    if (DFN[i]==LOW[i]) {
        Bcnt++;
        do {
            j=Stap[Stop--];
            instack[j]=false;
            Belong[j]=Bcnt;
        }
        while (j!=i);
    }
}
void solve() {
    int i;
    Stop=Bcnt=Dindex=0;
    memset(DFN,0,sizeof(DFN));
    for (i=1;i<=N;i++)
        if (!DFN[i])
            tarjan(i);
}
```

####KMP

####Astar

#数学 Mathematic
    数学最难之处在于现场推论，容易出现错误。
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
####取石子游戏
#####（一）巴什博奕（Bash Game）：

    只有一堆n个物品，两个人轮流从这堆物品中取物，规定每次至少取一个，最多取m个。最后取光者得胜。
    显然，如果n=m+1，那么由于一次最多只能取m个，所以，无论先取者拿走多少个，后取者都能够一次拿走剩余的物品，
    后者取胜。因此我们发现了如何取胜的法则：如果n=（m+1）r+s，（r为任意自然数，s≤m),那么先取者要拿走s个物品，
    如果后取者拿走k（≤m)个，那么先取者再拿走m+1-k个，结果剩下（m+1）（r-1）个，以后保持这样的取法，那么先取者肯定获胜。
    总之，要保持给对手留下（m+1）的倍数，就能最后获胜。

即，若n=k*(m+1)，则后取着胜，反之，存在先取者获胜的取法。

n%(m+1)==0. 先取者必败。
    这个游戏还可以有一种变相的玩法：两个人轮流报数，每次至少报一个，最多报十个，谁能报到100者胜。
    
    
从一堆100个石子中取石子，最后取完的胜。

#####（二）威佐夫博奕（Wythoff Game）：
    有两堆各若干个物品，两个人轮流从某一堆或同时从两堆中取同样多的物品，规定每次至少取一个，多者不限，最后取光者得胜。
    这种情况下是颇为复杂的。我们用（ak，bk）（ak ≤ bk ,k=0，1，2，...,n)表示两堆物品的数量并称其为局势，
    如果甲面对（0，0），那么甲已经输了，这种局势我们称为奇异局势。
    前几个奇异局势是：（0，0）、（1，2）、（3，5）、（4，7）、（6，10）、（8，13）、（9，15）、（11，18）、（12，20）。
    可以看出,a0=b0=0,ak是未在前面出现过的最小自然数,而 bk= ak + k，奇异局势有
如下三条性质：

    1。任何自然数都包含在一个且仅有一个奇异局势中。
    由于ak是未在前面出现过的最小自然数，所以有ak > ak-1 ，而 bk= ak + k > ak-1 + k-1 = bk-1 > ak-1 。所以性质1。成立。
    2。任意操作都可将奇异局势变为非奇异局势。
    事实上，若只改变奇异局势（ak，bk）的某一个分量，那么另一个分量不可能在其他奇异局势中，所以必然是非奇异局势。
    如果使（ak，bk）的两个分量同时减少，则由于其差不变，且不可能是其他奇异局势的差，因此也是非奇异局势。
    3。采用适当的方法，可以将非奇异局势变为奇异局势。
    

    假设面对的局势是（a,b），若 b = a，则同时从两堆中取走 a 个物体，就变为了奇异局势（0，0）
    如果a = ak ，b > bk，那么，取走b - bk个物体，即变为奇异局势
    如果 a = ak ， b < bk ,则同时从两堆中拿走 ak - ab - ak个物体,变为奇异局势（ ab - ak , ab - ak+ b - ak）
    如果a > ak ，b= ak + k,则从第一堆中拿走多余的数量a - ak 即可
    如果a < ak ，b= ak + k,分两种情况，第一种，a=aj （j < k）,从第二堆里面拿走 b - bj 即可
    第二种，a=bj （j < k）,从第二堆里面拿走 b - aj 即可。

从如上性质可知，两个人如果都采用正确操作，那么面对非奇异局势，先拿者必胜；反之，则后拿者取胜。

那么任给一个局势（a，b），怎样判断它是不是奇异局势呢？我们有如下公式：
    ak =[k（1+√5）/2]，bk= ak + k （k=0，1，2，...,n 方括号表示取整函数)
奇妙的是其中出现了黄金分割数（1+√5）/2 = 1。618...,因此,由ak，bk组成的矩形近似为黄金矩形，由于2/（1+√5）=（√5-1）/2，可以先求出j=[a（√5-1）/2]，若a=[j（1+√5）/2]，那么a = aj，bj = aj + j，若不等于，那么a = aj+1，bj+1 = aj+1+ j + 1，若都不是，那么就不是奇异局势。然后再按照上述法则进行，一定会遇到奇异局势。

#####（三）尼姆博奕（Nimm Game）：有三堆各若干个物品，两个人轮流从某一堆取任意多的物品，规定每次至少取一个，多者不限，最后取光者得胜。

    这种情况最有意思，它与二进制有密切关系，我们用（a，b，c）表示某种局势，
    首先（0，0，0）显然是奇异局势，无论谁面对奇异局势，都必然失败。
    第二种奇异局势是（0，n，n），只要与对手拿走一样多的物品，最后都将导致（0，0，0）。
    仔细分析一下，（1，2，3）也是奇异局势，无论对手如何拿，
        接下来都可以变为（0，n，n）的情形。

    计算机算法里面有一种叫做按位模2加，也叫做异或的运算，我们用符号（+）表示这种运算。
    这种运算和一般加法不同的一点是1+1=0。先看（1，2，3）的按位模2加的结果：

1 =二进制01
2 =二进制10
3 =二进制11 （+）
———————
0 =二进制00 （注意不进位）

    对于奇异局势（0，n，n）也一样，结果也是0。
    任何奇异局势（a，b，c）都有a（+）b（+）c =0。

如果我们面对的是一个非奇异局势（a，b，c），要如何变为奇异局势呢？假设 a < b< c,我们只要将 c 变为 a（+）b,即可,因为有如下的运算结果: a（+）b（+）(a（+）b)=(a（+）a)（+）(b（+）b)=0（+）0=0。要将c 变为a（+）b，只要从 c中减去 c-（a（+）b）即可。

    获胜情况对先取者进行讨论：
    异或结果为0，先取者必败，无获胜方法。后取者获胜；
    结果不为0，先取者有获胜的取法。

拓展： 任给N堆石子,两人轮流从任一堆中任取(每次只能取自一堆),取最后一颗石子的人获胜，问先取的人如何获胜？

根据上面所述，N个数异或即可。如果开始的时候T＝0，那么先取者必败，如果开始的时候T>0，那么只要每次取出石子使得T＝0，即先取者有获胜的方法。

 

 【综合一、三给出】

任给N堆石子,两人轮流从任一堆中任取(每次只能取自一堆),规定每方每次最多取K颗,取最后一颗石子的一方获胜.问先取的人如何获胜？

与上面的问题比，这个更复杂一些，我们可以这样做

    令Bi=Ai mod(K+1)
    定义T‘＝B1 xor B2 xor ... xor Bn
    如果T‘＝0 那么没有获胜可能，先取者必败
    如果T’>0 那么必然存在取的方法，使得T‘＝0，先取者有获胜的方法
    假设对方取了在Ai中取了r<=K个
    如果Ai中剩下的石子多于K 那么就在Ai中取走K+1-r个则Bi不变 T‘还是0
    如果Ai<=K 那么我们需要重新计算Bi和T‘ 按照上面的方法来做就可以了

####中国剩余定理
```	
推论1:
对于 a=ai  (mod ni) 的同余方程,有唯一解
下面说说由(a1, a2, ..., ak)求a的方法:
令 mi = n1*n2*...nk / ni;   ci = mi(mf  mod ni);   其中 mi*mf  mod ni = 1;
则 a = (a1*c1+a2*c2+...+ak*ck)      (mod n)      (注:由此等式可求a%n, 当n很大时)
剩余定理关键是mf的求法,如果理解了扩展欧几里得 ax+by=d, 就可以想到:
mi*mf  mod ni = 1 => mi*mf+ni*y=1;
```
```
int egcd(int a, int b, int &x, int &y) {
    int d;
    if (b == 0) {
        x = 1; y = 0; return a;
    } else {
        d = egcd(b, a % b, y, x);
        y -= a / b * x;
        return d;
    }
}
int lmes() {
    int i, tm=1, mf, y, ret=0, m;
    for (i = 0; i < nn; i++) tm *= n[i];
    for (i = 0; i < nn; i++) {
        m = tm / n[i];
        egcd(m, n[i], mf, y);
        ret += (a[i] * m * (mf % n[i])) % tm;
    }
    return (ret+tm) % tm;
}
```
####欧拉函数
```
若n是质数p的k次幂，φ(n)=p^k-p^(k-1)=(p-1)p^(k-1)
设n为正整数，以 φ(n)表示不超过n且与n互素的正整数的个数，
     称为n的欧拉函数值，这里函数φ：N→N，n→φ(n)称为欧拉函数。
欧拉函数是积性函数——若m,n互质，φ(mn)=φ(m)φ(n)。
特殊性质：当n为奇数时，φ(2n)=φ(n), 证明与上述类似。
```
####解析几何
	line: y-y0 = k(x-x0), k = tan(a), a = radius
	circle: (x-a)^2 + (y-b)^2 = r^2
	P(x0, y0), Ax+By+C=0, distance = abs(Ax0+By0+C)/sqrt(A^2+B^2)


####动态规划 方程
```c
剖分问题1-----石子合并
	f[i,j] = min(f[i,k]+f[k+1,j]+sum[i,j]);

LCS 最长公共子串
	f[i,j]={0                      (i=0)&(j=0);
    	   f[i-1,j-1]+1            (i>0,j>0,x=y[j]);
       	   max{f[i,j-1]+f[i-1,j]}} (i>0,j>0,x<>y[j]);

组合 递推
	C[I,j] = C[i-1,j]+C[I-1,j-1]
	C[I,0] = 1
最长公共子序列
	d[0,0]:=0;
	for i:=1 to n do
	  for j:=1 to m do
        if s1[i]=s2[j] then //等价于 if a[i,j]=1 then
          d[i,j]:=d[i-1,j-1]+1
        else begin
          if d[i-1,j]>d[i,j-1] then
            d[i,j]:=d[i-1,j]
          else d[i,j]:=d[i,j-1];
        end;

最长公共不下降子序列
	for(i=1;i<=l1;i++) {
  	    max=0;
  		for(j=1;j<=l2;j++)
    		if (b[j]<a[i]) {
      			max=MAX(max, f[j]);
    		} else if(b[j]==a[i]) {
      			f[j]=max+1;
      			ans=MAX(ans,f[j]); 
    		}
	}
	
01背包——> 简单
完全背包——> 简单
多重背包
```
