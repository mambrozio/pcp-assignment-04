#if !defined(graph_h)
#define graph_h

typedef int* Graph;

extern Graph graph;
extern int graph_size;

extern void graph_load(const char* path);
extern void graph_print(void);
extern void graph_free(void);

#define path(i, j) (graph[i * graph_size + j])

#endif
