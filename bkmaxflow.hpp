#ifndef BKMAXFLOW_HPP_INCLUDED
#define BKMAXFLOW_HPP_INCLUDED

#include<deque>
#include<vector>
#include<iostream>

namespace bk_maxflow
{
    struct edge;
    struct node;

    struct edge
    {
        node * head;
        edge * next;
        edge * sister;
        double r_cap;
        edge()
        {
            head = (node *)0;
            next = (edge *)0;
            sister = (edge *)0;
            r_cap = 0.0;
        }
	   ~edge()
	    {
	   		head = (node *)0;
	   		next = (edge *)0;
	   		sister = (edge *)0;
	    }
    };

    struct node
    {
        edge * first;
        edge * parent;
        node * next;
        int timestamp;
        int distance;
        bool in_sink;
        bool marked;
        bool in_changed_list;
        double r_cap;
        node()
        {
            first = (edge *)0;
            parent = (edge *)0;
            next = (node *)0;
            timestamp = 0;
            distance = 0;
            in_sink = false;
            marked = false;
            in_changed_list = false;
            r_cap = 0.0;
        }
       ~node()
        {
       		first = (edge *)0;
            parent = (edge *)0;
            next = (node *)0;
        }
    };


    enum terminal{FOREGROUND=0,BACKGROUND=1};
    static const double EPS = 1e-9;
    static const double INF = 1e100;
    class graph_cut
    {
        int num_nodes;
        int num_edges;

        double total_flow;

        int max_flow_iter;
        int time;

        edge to_terminal,to_orphan;

        node * active_first[2];
        node * active_last[2];
        std::deque<node *> orphans;
        node * nodes;
        edge * edges;

        void set_node_active(node * u)
        {
            if(!u->next)
            {
                if(active_last[1])
                    active_last[1]->next = u;
                else
                    active_first[1] = u;
                active_last[1] = u;
                u->next = u;
            }
        }

        node * get_next_active_node()
        {
            node * u = (node *)0;
            while(true)
            {
                u = active_first[0];
                if(!u)
                {
                    u = active_first[1];
                    active_first[0] = active_first[1];
                    active_last[0] = active_last[1];
                    active_first[1] = (node *)0;
                    active_last[1] = (node *)0;
                    if(!u)
                        return (node *)0;
                }
                if(u->next == u)
                {
                    active_first[0] = (node *)0;
                    active_last[0] = (node *)0;
                }
                else
                    active_first[0] = u->next;
                u->next = (node *)0;
                if(u->parent)
                    return u;
            }
            return (node *)0;
        }

        void add_orphan_at_front(node * u)
        {
            u->parent = &to_orphan;
            orphans.push_front(u);
        }

        void add_orphan_at_back(node * u)
        {
            u->parent = &to_orphan;
            orphans.push_back(u);
        }

        void max_flow_init()
        {
            active_first[0] = (node *)0;
            active_first[1] = (node *)0;
            active_last[0] = (node *)0;
            active_last[1] = (node *)0;
            orphans.clear();
            time = 0;
            for(int u=0;u<num_nodes;u++)
            {
                nodes[u].next = (node *)0;
                nodes[u].marked = false;
                nodes[u].in_changed_list = false;
                nodes[u].timestamp = time;
                if(nodes[u].r_cap > 0)
                {
                    nodes[u].in_sink = false;
                    nodes[u].parent = &to_terminal;
                    nodes[u].distance = 1;
                    set_node_active(&nodes[u]);
                }
                else if(nodes[u].r_cap < 0)
                {
                    nodes[u].in_sink = true;
                    nodes[u].parent = &to_terminal;
                    nodes[u].distance = 1;
                    set_node_active(&nodes[u]);
                }
                else nodes[u].parent = (edge *)0;
            }
        }

        void max_flow_reuse_trees_init()
        {
            node * node1 = (node *)0;
            node * node2 = (node *)0;
            node * queue_start = active_first[1];
            edge * current_edge = (edge *)0;

            active_first[0] = (node *)0;
            active_first[1] = (node *)0;
            active_last[0]  = (node *)0;
            active_last[1]  = (node *)0;

            orphans.clear();
            time++;

            while((node1 = queue_start))
            {
                queue_start = node1->next;
                if(queue_start == node1)
                    queue_start = (node *)0;
                node1->next = (node *)0;
                node1->marked = false;
                set_node_active(node1);
                if(node1->r_cap == 0)
                {
                    if(node1->parent)
                        add_orphan_at_back(node1);
                    continue;
                }
                if(node1->r_cap > 0)
                {
                    if(!node1->parent || node1->in_sink)
                    {
                        node1->in_sink = false;
                        for(current_edge = node1->first;
                            current_edge;current_edge = current_edge->next)
                        {
                            node2 = current_edge->head;
                            if(!node2->marked)
                            {
                                if(node2->parent == current_edge->sister)
                                    add_orphan_at_back(node2);
                                if(node2->parent && node2->in_sink && current_edge->r_cap > 0)
                                    set_node_active(node2);
                            }
                        }
                        node1->in_changed_list = true;
                    }
                }
                else
                {
                    if(!node1->parent || !node1->in_sink)
                    {
                        node1->in_sink = true;
                        for(current_edge = node1->first;
                            current_edge;current_edge = current_edge->next)
                        {
                            node2 = current_edge->head;
                            if(!node2->marked)
                            {
                                if(node2->parent == current_edge->sister)
                                    add_orphan_at_back(node2);
                                if(node2->parent && !node2->in_sink && current_edge->sister->r_cap > 0)
                                    set_node_active(node2);
                            }
                        }
                        node1->in_changed_list = true;
                    }
                }
                node1->parent = &to_terminal;
                node1->timestamp = time;
                node1->distance = 1;
            }

            while(orphans.size() > 0)
            {
                node * orphan = orphans.front();
                orphans.pop_front();
                if(orphan->in_sink)
                    process_sink_orphan(orphan);
                else
                    process_source_orphan(orphan);
            }

        }

        void augment(edge * middle)
        {
            node * current_node = (node *)0;
            edge * current_edge = (edge *)0;
            double delta = middle->r_cap;
            for(current_node = middle->sister->head;;current_node = current_edge->head)
            {
                current_edge = current_node->parent;
                if(current_edge == &to_terminal)
                    break;
                if(delta > current_edge->sister->r_cap)
                    delta = current_edge->sister->r_cap;
            }
            if(delta > current_node->r_cap)
                delta = current_node->r_cap;
            for(current_node = middle->head;;current_node = current_edge->head)
            {
                current_edge = current_node->parent;
                if(current_edge == &to_terminal)
                    break;
                if(delta > current_edge->r_cap)
                    delta = current_edge->r_cap;
            }
            if(delta > -current_node->r_cap)
                delta = -current_node->r_cap;

            middle->sister->r_cap += delta;
            middle->r_cap -=delta;
            for(current_node = middle->sister->head;;current_node = current_edge->head)
            {
                current_edge = current_node->parent;
                if(current_edge == &to_terminal)
                    break;
                current_edge->r_cap += delta;
                current_edge->sister->r_cap -= delta;
                if(current_edge->sister->r_cap == 0)
                    add_orphan_at_front(current_node);
            }
            current_node->r_cap -= delta;
            if(current_node->r_cap == 0)
                add_orphan_at_front(current_node);

            for(current_node = middle->head;;current_node = current_edge->head)
            {
                current_edge = current_node->parent;
                if(current_edge == &to_terminal)
                    break;
                current_edge->r_cap -= delta;
                current_edge->sister->r_cap += delta;
                if(current_edge->r_cap == 0)
                    add_orphan_at_front(current_node);
            }
            current_node->r_cap += delta;
            if(current_node->r_cap == 0)
                add_orphan_at_front(current_node);

            total_flow += delta;
        }

        void process_source_orphan(node * orphan)
        {
            edge * best_edge = (edge *) 0;
            int min_distance = -1;
            for(edge * orphan_edge = orphan->first; orphan_edge; orphan_edge = orphan_edge->next)
            if(orphan_edge->sister->r_cap != 0)
            {
                node * node_aux = orphan_edge->head;
                edge * parent_edge = node_aux->parent;
                if(parent_edge && !node_aux->in_sink)
                {
                    int distance = 0;
                    while(true)
                    {
                        if(node_aux->timestamp == time)
                        {
                            distance += node_aux->distance;
                            break;
                        }
                        parent_edge = node_aux->parent;
                        distance++;
                        if(parent_edge == &to_terminal)
                        {
                            node_aux->timestamp = time;
                            node_aux->distance = 1;
                            break;
                        }
                        if(parent_edge == &to_orphan)
                        {
                            distance = -1;
                            break;
                        }
                        node_aux = parent_edge->head;
                    }
                    if(distance >= 0)
                    {
                        if(distance < min_distance || min_distance < 0)
                        {
                            best_edge = orphan_edge;
                            min_distance = distance;
                        }
                        for(node_aux = orphan_edge->head;node_aux->timestamp != time; node_aux = node_aux->parent->head)
                        {
                            node_aux->timestamp = time;
                            node_aux->distance = distance;
                            distance--;
                        }
                    }
                }
            }
            orphan->parent = best_edge;
            if(best_edge)
            {
                orphan->timestamp = time;
                orphan->distance = min_distance + 1;
            }
            else
            {
                orphan->in_changed_list = true;
                for(edge * orphan_edge = orphan->first;orphan_edge;orphan_edge = orphan_edge->next)
                {
                    node * node_aux = orphan_edge->head;
                    edge * parent_edge = node_aux->parent;
                    if(!node_aux->in_sink && parent_edge)
                    {
                        if(orphan_edge->sister->r_cap != 0)
                            set_node_active(node_aux);
                        if(parent_edge != &to_terminal && parent_edge != &to_orphan
                        && parent_edge->head == orphan)
                            add_orphan_at_back(node_aux);
                    }
                }
            }
        }

        void process_sink_orphan(node * orphan)
        {
            edge * best_edge = (edge *)0;
            int min_distance = -1;
            for(edge * orphan_edge = orphan->first;orphan_edge;orphan_edge = orphan_edge->next)
            if(orphan_edge->r_cap != 0)
            {
                node * node_aux = orphan_edge->head;
                edge * parent_edge = node_aux->parent;
                if(node_aux->in_sink && parent_edge)
                {
                    int distance = 0;
                    while(true)
                    {
                        if(node_aux->timestamp == time)
                        {
                            distance += node_aux->distance;
                            break;
                        }
                        parent_edge = node_aux->parent;
                        distance++;
                        if(parent_edge == &to_terminal)
                        {
                            node_aux->timestamp = time;
                            node_aux->distance = 1;
                            break;
                        }
                        if(parent_edge == &to_orphan)
                        {
                            distance = -1;
                            break;
                        }
                        node_aux = parent_edge->head;
                    }
                    if(distance >= 0)
                    {
                        if(distance < min_distance || min_distance == -1)
                        {
                            best_edge = orphan_edge;
                            min_distance = distance;
                        }
                        for(node_aux = orphan_edge->head;node_aux->timestamp != time;node_aux = node_aux->parent->head)
                        {
                            node_aux->timestamp = time;
                            node_aux->distance = distance;
                            distance--;
                        }
                    }
                }
            }
            orphan->parent = best_edge;
            if(best_edge)
            {
                orphan->timestamp = time;
                orphan->distance = min_distance + 1;
            }
            else
            {
                orphan->in_changed_list = true;
                for(edge * orphan_edge = orphan->first;orphan_edge;orphan_edge = orphan_edge->next)
                {
                    node * node_aux = orphan_edge->head;
                    edge * parent_edge = node_aux->parent;
                    if(node_aux->in_sink && parent_edge)
                    {
                        if(orphan_edge->r_cap != 0)
                            set_node_active(node_aux);
                        if(parent_edge != &to_terminal && parent_edge != &to_orphan
                        && parent_edge->head == orphan)
                            add_orphan_at_back(node_aux);
                    }
                }
            }
        }

        public:
            graph_cut()
            {
                num_nodes = 0;
                num_edges = 0;
                max_flow_iter = 0;
                time = 0;
                total_flow = 0.0;
            }
		   ~graph_cut()
		    {
		   	delete [] edges;
		   	delete [] nodes;
		   	orphans.clear();
		   	active_first[0] = (node *)0;
		   	active_first[1] = (node *)0;
		   	active_last[0] = (node *)0;
		   	active_last[0] = (node *)0;
			}
            void set_terminal_weights(int u_idx,double to_source,double to_sink)
            {
                double delta = nodes[u_idx].r_cap;
                if(delta > 0)
                    to_source += delta;
                else
                    to_sink   -= delta;
                total_flow += (to_source < to_sink) ? to_source : to_sink;
                nodes[u_idx].r_cap = to_source - to_sink;
            }

            void set_edge_weight(int u_idx,int v_idx,double weight_uv,double weight_vu)
            {
                //edges.push_back(edge());
                num_edges += 2;
                //edges.push_back(edge());
                //num_edges++;

                edges[num_edges-2].sister = &edges[num_edges-1];
                edges[num_edges-1].sister = &edges[num_edges-2];

                edges[num_edges-2].next = nodes[u_idx].first;
                nodes[u_idx].first = &edges[num_edges-2];

                edges[num_edges-1].next = nodes[v_idx].first;
                nodes[v_idx].first = &edges[num_edges-1];

                edges[num_edges-2].head = &nodes[v_idx];
                edges[num_edges-1].head = &nodes[u_idx];

                edges[num_edges-2].r_cap = weight_uv;
                edges[num_edges-1].r_cap = weight_vu;
            }
            double get_maxflow(bool reuse_trees,std::vector<int> & changed_nodes)
            {
                if(max_flow_iter == 0)
                    reuse_trees = false;
                if(reuse_trees)
                    max_flow_reuse_trees_init();
                else
                    max_flow_init();

                node * current_node = (node *)0;
                edge * current_edge = (edge *)0;

                while(true)
                {
                    node * active_node = current_node;
                    if(active_node)
                    {
                        active_node->next = (node *)0;
                        if(!active_node->parent)
                            active_node = (node *)0;
                    }
                    else
                    {
                        active_node = get_next_active_node();

                    }
                    if(!active_node)
                            break;
                    if(!active_node->in_sink)
                    {
                        for(current_edge = active_node->first;
                            current_edge; current_edge = current_edge->next)
                            if(current_edge->r_cap > EPS)
                            {
                                node * head_node = current_edge->head;
                                if(!head_node->parent)
                                {
                                    head_node->in_sink = false;
                                    head_node->parent = current_edge->sister;
                                    head_node->timestamp = active_node->timestamp;
                                    head_node->distance = active_node->distance + 1;
                                    set_node_active(head_node);
                                    head_node->in_changed_list = true;
                                }
                                else if(head_node->in_sink)
                                    break;
                                else if(head_node->timestamp <= active_node->timestamp
                                    && head_node->distance > active_node->distance)
                                {
                                    head_node->parent = current_edge->sister;
                                    head_node->timestamp = active_node->timestamp;
                                    head_node->distance = active_node->distance + 1;
                                }
                            }
                    }
                    else
                    {
                        for(current_edge = active_node->first;
                        current_edge;current_edge = current_edge->next)
                            if(current_edge->sister->r_cap > EPS)
                            {
                                node * head_node = current_edge->head;
                                if(!head_node->parent)
                                {
                                    head_node->in_sink = true;
                                    head_node->parent = current_edge->sister;
                                    head_node->timestamp = active_node->timestamp;
                                    head_node->distance = active_node->distance + 1;
                                    set_node_active(head_node);
                                    head_node->in_changed_list = true;
                                }
                                else if(!head_node->in_sink)
                                {
                                    current_edge = current_edge->sister;
                                    break;
                                }
                                else if(head_node->timestamp <= active_node->timestamp
                                    && head_node->distance > active_node->distance)
                                {
                                    head_node->parent = current_edge->sister;
                                    head_node->timestamp = active_node->timestamp;
                                    head_node->distance = active_node->distance + 1;
                                }
                            }
                    }

                    time++;

                    if(current_edge)
                    {
                        active_node->next = active_node;
                        current_node = active_node;

                        augment(current_edge);

                        while(orphans.size() > 0)
                        {
                            node * orphan = orphans.front();
                            orphans.pop_front();
                            if(orphan->in_sink)
                                process_sink_orphan(orphan);
                            else
                                process_source_orphan(orphan);
                        }
                    }
                    else current_node = (node *)0;
                }
                max_flow_iter++;
                if(changed_nodes.size()>0)
                {
                    changed_nodes.clear();
                    for(int i=0;i<num_nodes;i++)
                        if(nodes[i].in_changed_list)
                            changed_nodes.push_back(i);
                }
                return total_flow;
            }

            terminal get_terminal(int u_idx)
            {
                if(nodes[u_idx].parent)
                    return nodes[u_idx].in_sink ? BACKGROUND:FOREGROUND;
                else return BACKGROUND;
            }

            int get_num_edges()
            {
                return this->num_edges;
            }

            int get_num_nodes()
            {
                return this->num_nodes;
            }

            bool connected(int u_idx,int v_idx)
            {
                //std::cout<<"checking"<<u_idx<<" "<<v_idx<<std::endl;
                for(edge * current_edge = nodes[u_idx].first; current_edge ;current_edge = current_edge->next)
            		{
            		//W.printf("%x",current_edge->head);
            		if(current_edge->head == &nodes[v_idx])return true;
            		}
                return false;
            }
            void add_nodes(int num)
            {
            	num_nodes += num;
                nodes = new node [num_nodes];
            }

            void add_edges(int num)
            {
                edges = new edge [num];
            }

            void mark_node(int u_idx)
            {
                if(!nodes[u_idx].next)
                {
                    if(active_last[1])
                        active_last[1]->next = &nodes[u_idx];
                    else
                        active_first[1] = &nodes[u_idx];
                    active_last[1] = &nodes[u_idx];
                    nodes[u_idx].next = &nodes[u_idx];
                }
                nodes[u_idx].marked = true;
            }
    };
};

#endif // BKMAXFLOW_HPP_INCLUDED
