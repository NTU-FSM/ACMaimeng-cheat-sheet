ACM-cheat-sheet
===============

#基础 Basic
####Buffered Input
```java
BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
    int T = Integer.parseInt(br.readLine());
    String[] data;
    while (T-- > 0) {
        data = br.readLine().split(" ");
        n = Integer.parseInt(data[0]);
        m = Integer.parseInt(data[1]);
```

####String (Regular Expression)
```java
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

####Set(?)
```java
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
```

####RMQ(Range Minimum Query)
[for LCA(Lowest Common Ancestor)](http://community.topcoder.com/tc?module=Static&d1=tutorials&d2=lowestCommonAncestor)


#算法 algorithm
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
```java
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
