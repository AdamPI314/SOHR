///////////////////////////////////////////////////////////////////////////////////////////////////
//This is a boost graph template with bundled property
///////////////////////////////////////////////////////////////////////////////////////////////////
/************************************************************************/
//http://stackoverflow.com/questions/671714/modifying-vertex-properties-in-a-boostgraph
/*   Examples:
struct VertexProperties_graph {
int i;
};

struct EdgeProperties_graph {
};

typedef Graph_bundled<VertexProperties_graph, EdgeProperties_graph> MyGraph;

MyGraph g;

VertexProperties_graph vp;
vp.i = 42;

MyGraph::Vertex v = g.AddVertex(vp);

g.properties(v).i = 23;*/
/************************************************************************/

#ifndef GRAPH_BUNDLED_H_
#define GRAPH_BUNDLED_H_

#include <boost/config.hpp>
#include <boost/version.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/static_assert.hpp>
#include <algorithm> //for std::find_if

using namespace boost;

/* definition of basic boost::graph properties */
enum vertex_properties_t { vertex_properties };
enum edge_properties_t { edge_properties };
//enum graph_properties_t { graph_properties };

namespace boost {
	BOOST_INSTALL_PROPERTY(vertex, properties);
	BOOST_INSTALL_PROPERTY(edge, properties);
	//BOOST_INSTALL_PROPERTY(graph, properties);
}


/* the graph base class template */
template < typename VERTEXPROPERTIES, typename EDGEPROPERTIES >
class Graph_bundled
{
public:
	/* an adjacency_list as we need it */
	/*
	An adjacency_list doesn't have an edge index associated with it, only a vertex index.
	Which is quite logical once you think about how the graph is stored.
	To have an edge index, you need to manually add it to the graph description, and then manually handle it.
	*/
	typedef adjacency_list<
		vecS, // allow parallel edges
		//multisetS, //Specify boost::multisetS as the OutEdgeListS template argument to adjacency_list, this will enable an extra function
		//called edge_range which returns a range of iterators for all the "parallel" edges coming out of u and going into v
		//listS, // doesn't allow parallel edges
		vecS, // vertex container
		bidirectionalS, // directed graph
		/*property<vertex_properties_t, VERTEXPROPERTIES>,*/
		//property<edge_properties_t, EDGEPROPERTIES>

		property<vertex_index_t, std::size_t, property<vertex_properties_t, VERTEXPROPERTIES> >,
		property<edge_index_t, std::size_t, property<edge_properties_t, EDGEPROPERTIES> > //multiple properties
		//edge_index_t is not installed for default, if you want to use get() method to get edge_index, you have to install it explicitly
	> GraphContainer;

	/* a bunch of graph-specific typedefs */
	typedef typename graph_traits<GraphContainer>::vertex_descriptor Vertex;
	typedef std::pair<Vertex, Vertex> VertexPair;
	typedef typename graph_traits<GraphContainer>::edge_descriptor Edge;
	typedef std::pair<Edge, Edge> EdgePair;

	typedef typename graph_traits<GraphContainer>::vertex_iterator vertex_iter;
	typedef typename graph_traits<GraphContainer>::edge_iterator edge_iter;
	typedef typename graph_traits<GraphContainer>::adjacency_iterator adjacency_iter;
	typedef typename graph_traits<GraphContainer>::out_edge_iterator out_edge_iter;

	typedef typename graph_traits<GraphContainer>::degree_size_type degree_t;

	typedef std::pair<adjacency_iter, adjacency_iter> adjacency_vertex_range_t;
	typedef std::pair<out_edge_iter, out_edge_iter> out_edge_range_t;
	typedef std::pair<vertex_iter, vertex_iter> vertex_range_t;
	typedef std::pair<edge_iter, edge_iter> edge_range_t;

	//typedef typename property_map<GraphContainer, edge_index_t>::type edge_index_map_t;

	//typedef for reversed graph
	typedef reverse_graph<GraphContainer> reverse_graph_t;
	typedef typename graph_traits<reverse_graph_t>::vertex_descriptor reverse_vertex_t;
	typedef typename graph_traits<reverse_graph_t>::edge_descriptor reverse_edge_t;

public:
	//visitor for shortest path algotithm, both dijkstra and bellman_ford
	//here we just implement the "edge_relaxed" function so that it could record the edge on the shortest path
	//class my_visitor : public default_bfs_visitor
	//class my_visitor : public dijkstra_visitor<null_visitor>
	class my_visitor : public base_visitor<null_visitor>
	{
	protected:
		//predecessors vector
		std::vector<Vertex>& predecessors;
		//coresponding edges vector of predecessors, edge from current vertex<--predessor coressponds a edge
		std::vector<Edge>& edge_v;

	public:
		template <typename Graph>
		my_visitor(Graph g, std::vector<Vertex>& p, std::vector<Edge>& e_v) :predecessors(p), edge_v(e_v)
		{
			predecessors.resize(num_vertices(g));
			edge_v.resize(num_vertices(g));
		}
	public:
		template <typename Vertex, typename Graph>
		void initialize_vertex(const Vertex &s, const Graph &g) const {}
		template <typename Vertex, typename Graph>
		void discover_vertex(const Vertex &s, const Graph &g) const {}
		template <typename Vertex, typename Graph>
		void examine_vertex(const Vertex &s, const Graph &g) const {}
		template <typename Edge, typename Graph>
		void examine_edge(const Edge &e, const Graph &g) const {}
		//no const statement for this function, gonna change the class itself
		template <typename Edge, typename Graph>
		void edge_relaxed(const Edge &e, const Graph &g)
		{
			auto s = source(e, g);
			auto t = target(e, g);
			//my own way to create predecessor map, or says predecessor vector
			predecessors[t] = s;
			edge_v[t] = e;
			//std::cout << s << "-->" << t << std::endl;

		}
		template <typename Edge, typename Graph>
		void edge_not_relaxed(const Edge &e, const Graph &g) const {}
		template <typename Vertex, typename Graph>
		void finish_vertex(const Vertex &s, const Graph &g) const {}
		//bellman_ford needs this two extra functions
		template <typename Edge, typename Graph>
		void edge_minimized(const Edge &e, const Graph &g) const {}
		template <typename Edge, typename Graph>
		void edge_not_minimized(const Edge &e, const Graph &g) const {}


	};
public:
	//visitor for shortest path algotithm, both dijkstra and bellman_ford
	//here we just implement the "edge_relaxed" function so that it could record the edge on the shortest path
	//class my_visitor : public default_bfs_visitor
	//class my_visitor : public dijkstra_visitor<null_visitor>
	//for reverse graph
	class my_visitor_r : public base_visitor<null_visitor>
	{
	protected:
		//predecessors vector
		std::vector<reverse_vertex_t>& predecessors;
		//coresponding reverse_edge_ts vector of predecessors, reverse_edge_t from current reverse_vertex_t<--predessor coressponds a reverse_edge_t
		std::vector<reverse_edge_t>& reverse_edge_v;

	public:
		template <typename reverse_graph_t>
		my_visitor_r(reverse_graph_t g, std::vector<reverse_vertex_t>& p, std::vector<reverse_edge_t>& e_v) :predecessors(p), reverse_edge_v(e_v)
		{
			predecessors.resize(num_vertices(g));
			reverse_edge_v.resize(num_vertices(g));
		}
	public:
		template <typename reverse_vertex_t, typename reverse_graph_t>
		void initialize_vertex(const reverse_vertex_t &s, const reverse_graph_t &g) const {}
		template <typename reverse_vertex_t, typename reverse_graph_t>
		void discover_vertex(const reverse_vertex_t &s, const reverse_graph_t &g) const {}
		template <typename reverse_vertex_t, typename reverse_graph_t>
		void examine_vertex(const reverse_vertex_t &s, const reverse_graph_t &g) const {}
		template <typename reverse_edge_t, typename reverse_graph_t>
		void examine_edge(const reverse_edge_t &e, const reverse_graph_t &g) const {}
		//no const statement for this function, gonna change the class itself
		template <typename reverse_edge_t, typename reverse_graph_t>
		void edge_relaxed(const reverse_edge_t &e, const reverse_graph_t &g)
		{
			auto s = source(e, g);
			auto t = target(e, g);
			//my own way to create predecessor map, or says predecessor vector
			predecessors[t] = s;
			reverse_edge_v[t] = e;

		}
		template <typename reverse_edge_t, typename reverse_graph_t>
		void edge_not_relaxed(const reverse_edge_t &e, const reverse_graph_t &g) const {}
		template <typename reverse_vertex_t, typename reverse_graph_t>
		void finish_vertex(const reverse_vertex_t &s, const reverse_graph_t &g) const {}
		//bellman_ford needs this two extra functions
		template <typename reverse_edge_t, typename reverse_graph_t>
		void edge_minimized(const reverse_edge_t &e, const reverse_graph_t &g) const {}
		template <typename reverse_edge_t, typename reverse_graph_t>
		void edge_not_minimized(const reverse_edge_t &e, const reverse_graph_t &g) const {}

	};

	/* constructors etc. */
	Graph_bundled() {}

	Graph_bundled(const Graph_bundled& g) : graph(g.graph) {}

	virtual ~Graph_bundled() {}

	/* structure modification methods */
	void Clear()
	{
		graph.clear();
	}

	Vertex AddVertex(const VERTEXPROPERTIES& prop)
	{
		Vertex v = add_vertex(graph);
		properties(v) = prop;
		return v;
	}

	void RemoveVertex(const Vertex& v)
	{
		clear_vertex(v, graph);
		remove_vertex(v, graph);
	}

	Edge AddEdge(const Vertex& v1, const Vertex& v2, const EDGEPROPERTIES& prop_12)
	{
		/* TODO: maybe one wants to check if this edge could be inserted */
		Edge addedEdge = add_edge(v1, v2, graph).first;

		properties(addedEdge) = prop_12;

		return addedEdge;
	}

	EdgePair AddEdge(const Vertex& v1, const Vertex& v2, const EDGEPROPERTIES& prop_12, const EDGEPROPERTIES& prop_21)
	{
		/* TODO: maybe one wants to check if this edge could be inserted */
		Edge addedEdge1 = add_edge(v1, v2, graph).first;
		Edge addedEdge2 = add_edge(v2, v1, graph).first;

		properties(addedEdge1) = prop_12;
		properties(addedEdge2) = prop_21;

		return EdgePair(addedEdge1, addedEdge2);
	}

	/* property access */
	VERTEXPROPERTIES& properties(const Vertex& v)
	{
		typename property_map<GraphContainer, vertex_properties_t>::type param = get(vertex_properties, graph);
		return param[v];
	}

	const VERTEXPROPERTIES& properties(const Vertex& v) const
	{
		typename property_map<GraphContainer, vertex_properties_t>::const_type param = get(vertex_properties, graph);
		return param[v];
	}

	EDGEPROPERTIES& properties(const Edge& v)
	{
		typename property_map<GraphContainer, edge_properties_t>::type param = get(edge_properties, graph);
		return param[v];
	}

	const EDGEPROPERTIES& properties(const Edge& v) const
	{
		typename property_map<GraphContainer, edge_properties_t>::const_type param = get(edge_properties, graph);
		return param[v];
	}

	const EDGEPROPERTIES& properties(const reverse_edge_t& v) const
	{
		//reversed graph
		reverse_graph_t reversed_graph(this->graph);
		typename property_map<reverse_graph_t, edge_properties_t>::const_type param = get(edge_properties, reversed_graph);
		return param[v];
	}

	/* selectors and properties */
	const GraphContainer& getGraph() const
	{
		return graph;
	}

	vertex_range_t getVertices() const
	{
		return vertices(graph);
	}

	std::size_t get_num_vertices() const {
		vertex_range_t vp;
		int count = 0;
		for (vp = vertices(graph); vp.first != vp.second; ++vp.first)
			count++;
		return count;
	}

	edge_range_t getEdges() const
	{
		return edges(graph);
	}

	std::size_t getEdgeCount() const
	{
		return num_edges(graph);
	}

	std::size_t get_num_edges() const {
		edge_range_t ep;
		int count = 0;
		for (ep = edges(graph); ep.first != ep.second; ++ep.first)
			count++;
		return count;
	}

	adjacency_vertex_range_t getAdjacentVertices(const Vertex& v) const
	{
		return adjacent_vertices(v, graph);
	}

	std::size_t getInDegree(const Vertex&v) const {
		return in_degree(v, graph);
	}

	std::size_t getOutDegree(const Vertex&v) const {
		return out_degree(v, graph);
	}

	out_edge_range_t getOutEdges(const Vertex&v) const {
		return out_edges(v, graph);
	}

	std::size_t getVertexCount() const
	{
		return num_vertices(graph);
	}

	std::size_t getVertexDegree(const Vertex& v) const
	{
		return out_degree(v, graph);
	}

	/* operators */
	Graph_bundled& operator=(const Graph_bundled &rhs)
	{
		graph = rhs.graph;
		return *this;
	}

	//find the first edge with specific edge_index
	edge_iter find_edge_with_an_edge_index(size_t x) {
		edge_iter iter_beg, iter_end;
		boost::tie(iter_beg, iter_end) = getEdges();

		//return std::find_if(iter_beg, iter_end, check_x(x));		
		for (; iter_beg != iter_end; ++iter_beg) {
			if (properties(*iter_beg).edge_index == x) {
				break;
			}
		}
		return iter_beg;
	}

protected:
	GraphContainer graph;

};


#endif
