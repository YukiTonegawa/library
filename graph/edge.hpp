#ifndef _EDGE_H_
#define _EDGE_H_
/*
template<typename edge_weight>
struct edge_base{
  using weight = edge_weight;
  int to();
  int from();
  int id();
  weight wei();
  static weight z();
  edge_base<weight> reverse();
};
*/

template<typename edge_weight>
struct simple_edge{
  using weight = edge_weight;
  int s, t;
  simple_edge(): s(-1), t(-1){}
  simple_edge(int a, int b): s(a), t(b){}
  int to(){return t;}
  int from(){return s;}
  int id(){return -1;}
  weight wei(){return 1;}
  static weight z(){return 0;}
  simple_edge<weight> reverse(){return simple_edge<weight>{t, s};}
};

template<typename edge_weight>
struct weighted_edge{
  using weight = edge_weight;
  int s, t;
  weight w;
  weighted_edge(): s(-1), t(-1), w(0){}
  weighted_edge(int a, int b, weight c): s(a), t(b), w(c){}
  int to(){return t;}
  int from(){return s;}
  int id(){return -1;}
  weight wei(){return w;}
  static weight z(){return 0;}
  weighted_edge<weight> reverse(){return weighted_edge<weight>{t, s, w};}
};

template<typename edge_weight>
struct labeled_edge{
  using weight = edge_weight;
  int s, t, i;
  labeled_edge(): s(-1), t(-1), i(-1){}
  labeled_edge(int a, int b, int i): s(a), t(b), i(i){}
  int to(){return t;}
  int from(){return s;}
  int id(){return i;}
  weight wei(){return 1;}
  static weight z(){return 0;}
  labeled_edge<weight> reverse(){return labeled_edge<weight>{t, s, i};}
};

template<typename edge_weight>
struct weighted_labeled_edge{
  using weight = edge_weight;
  int s, t;
  weight w;
  int i;
  weighted_labeled_edge(): s(-1), t(-1), w(0), i(-1){}
  weighted_labeled_edge(int a, int b, weight w, int i): s(a), t(b), w(w), i(i){}
  int to(){return t;}
  int from(){return s;}
  int id(){return i;}
  weight wei(){return w;}
  static weight z(){return 0;}
  weighted_labeled_edge<weight> reverse(){return weighted_labeled_edge<weight>{t, s, w, i};}
};
#endif