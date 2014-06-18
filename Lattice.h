#ifndef _LATTICE_H
#define _LATTICE_H

#include <vector>
#include <map>
#include <string>
#include <iostream>

class Edge;

class Vertex {
    public:
        size_t idx;
        std::vector<size_t> uc; //!< Unit cell coordinates
        std::vector<double> pos; //!< Position within unit cell
        std::vector<Edge*> edges;
        Vertex(size_t idx, const std::vector<size_t>& uc, const std::vector<double>& pos)
            : idx(idx), uc(uc), pos(pos) 
        {}
        ~Vertex()
        {}
};

class Edge {
    public:
        Vertex* first;
        Vertex* second;
        double Jprop;
        Edge(Vertex* first, Vertex* second, double Jprop=1)
            : first(first), second(second), Jprop(Jprop)
        {
            first->edges.push_back(this);
            second->edges.push_back(this);
        }
        const Vertex* GetOther(const Vertex* v) const
        {
            return v==first? second : first;
        }
        ~Edge()
        {}
};

class Lattice {
    protected:
        std::vector<const Vertex*> vertices;
        std::vector<const Edge*> edges;
        size_t m_Lx; //!< number of unit cells in x direction
        size_t m_Ly; //!< number of unit cells in y direction
    public:
        Lattice(size_t Lx,size_t Ly)
            :m_Lx(Lx), m_Ly(Ly) {}
        ~Lattice()
        {
            for(size_t e=0;e<edges.size();++e)
                delete edges[e];
            for(size_t v=0;v<vertices.size();++v)
                delete vertices[v];
        }
        size_t GetNv() const
        {
            return vertices.size();
        }
        size_t GetLx() const
        {
            return m_Lx;
        }
        size_t GetLy() const
        {
            return m_Ly;
        }
        const std::vector<const Vertex*>& GetVertices() const
        {
            return vertices;
        }
        const std::vector<const Edge*>& GetEdges() const
        {
            return edges;
        }
        virtual std::string str(std::vector<std::string> st) const =0;
        friend std::ostream& operator<<(std::ostream& out,const Lattice& lat)
        {
            out<<(&lat)->str(std::vector<std::string>(lat.vertices.size()));
            return out;
        }
};

#endif//LATTICE_H
