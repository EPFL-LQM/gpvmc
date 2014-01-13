#ifndef _LATTICE_H
#define _LATTICE_H

#include <vector>
#include <map>
#include <string>

class Edge;

class Vertex {
    public:
        std::vector<double> pos;
        std::vector<Edge*> edges;
        Vertex(std::vector<double> pos, size_t state=0)
            : pos(pos)
        {}
        ~Vertex()
        {}
};

class Edge {
    public:
        Vertex* first;
        Vertex* second;
        std::map<std::string,double> prop;
        Edge(Vertex* first, Vertex* second, std::string propname="", double prop=0)
            : first(first), second(second)
        {
            if(propname.size())
                this->prop[propname]=prop;
            first->edges.push_back(this);
            second->edges.push_back(this);
        }
        ~Edge()
        {}
};

class Lattice {
    public:
        std::vector<Vertex*> vertices;
        std::vector<Edge*> edges;
        Lattice()
        {}
        ~Lattice()
        {
            for(size_t e=0;e<edges.size();++e)
                delete edges[e];
            for(size_t v=0;v<vertices.size();++v)
                delete vertices[v];
        }
};

#endif//LATTICE_H
