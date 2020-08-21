var EulerHeuristic = (function () {
  function EulerHeuristic(graph, graphA) {
    this._graph = graph;
    this._graphA = graphA;
    this._tour = [];
  }

  EulerHeuristic.prototype.getTour = function () {
    return this._tour;
  };
  EulerHeuristic.prototype.containsArc = function (id) {
    for (var ind in this._tour) {
      if (this._tour[ind].getId() == id) return 1;
    }
    return 0;
  };

  EulerHeuristic.prototype.run = function () {
    var tabEueler = this._graphA.cycleEuler();
    var T_arc = this._graph.getTarc();
    var tour = [];
    // initialisation des deux colonnes pères et fils de T_arc
    for (var i = 0; i < T_arc.length; i++) {
      T_arc[i][2] = new Array();
      T_arc[i][3] = new Array();
    }
    var arc;
    // parcours dans le sens direct
    for (var ind in tabEueler) {
      arc = this._graphA.getEdges()[tabEueler[ind]];
      T_arc = this._graphA.setPeresFils(tour, arc, T_arc);
    }
    // parcours dans le sens inverse
    for (var i = tabEueler.length - 1; i >= 0; i--) {
      arc = this._graphA.getEdges()[tabEueler[i]];
      arc = this._graphA.findEdge(arc.getDest(), arc.getSrc());
      T_arc = this._graphA.setPeresFils(tour, arc, T_arc);
    }
    // affectation de la solution
    this._tour = tour;
  };

  return EulerHeuristic;
})();

var Tour = (function () {
  function Tour(graphA) {
    this._graphA = graphA;
    this._tour = [];
  }

  Tour.prototype.size = function () {
    return this._tour.length;
  };

  Tour.prototype.contains = function (noeudA) {
    for (var tourIndex in this._tour) {
      if (noeudA.isEqual(this._tour[tourIndex])) {
        return true;
      }
    }

    return false;
  };

  Tour.prototype.addNoeudA = function (noeudA) {
    this._tour.push(noeudA);
  };

  Tour.prototype.getNoeudA = function (tourIndex) {
    return this._tour[tourIndex];
  };

  return Tour;
})();

var Eval = (function () {
  function Eval(graph, solEuler, niveaux) {
    this._graph = graph;
    this._solEuler = solEuler;
    this._tabRoutage = [];
    this._niveaux = niveaux;
  }

  Eval.prototype.addACOsol = function (solACO) {
    this._solACO = solACO;
  };

  // Eval.prototype.cout=function(srce,dest){
  // var cmin=0;
  // var c1=0;
  // var c2=0;
  // var s,d;

  // for(var i=0; i< this._graph.getDimension();i++){
  // s=Math.min(srce.getCoord(i),dest.getCoord(i));
  // d=Math.max(srce.getCoord(i),dest.getCoord(i));
  // c1=Math.abs(srce.getCoord(i)-dest.getCoord(i));

  // c2=(s%(this._graph.getBase()-1))+1+Math.abs((this._graph.getBase()-1)-d);
  // cmin=cmin+Math.min(c1,c2);
  // }
  // return cmin;
  // };

  // on calcul cette longeur en utilisant la hamming distance ou le xor bit a bit

  Eval.prototype.cout = function (srce, dest) {
    var cmin = 0;

    for (var i = 0; i < this._graph.getDimension(); i++) {
      cmin = cmin + (srce.getCoord(i) ^ dest.getCoord(i));
    }
    return cmin;
  };

  Eval.prototype.calculTableRoutage = function () {
    var srce;
    var dest;
    var v;
    var index;
    var nbNoeuds = this._graph.getNoeuds().length;
    //initialisation de la table de routage
    for (var i = 0; i < nbNoeuds; i++) {
      this._tabRoutage[i] = new Array(nbNoeuds);
    }
    for (var i = 0; i < nbNoeuds; i++) {
      for (var j = 0; j < nbNoeuds; j++) {
        this._tabRoutage[i][j] = new Array();
      }
    }
    //remplissage de la table de routage avec les dépendances permises par la méthode utilisée
    for (var i = 0; i < nbNoeuds; i++) {
      srce = this._graph.getNoeuds()[i];
      for (var j = 0; j < nbNoeuds; j++) {
        dest = this._graph.getNoeuds()[j];
        if (
          srce.getId() != dest.getId() &&
          this._graph.estVoisin(srce, dest) == 0
        ) {
          for (var d = 0; d < this._graph.getTableVoisins()[i].length; d++) {
            v = this._graph.getTableVoisins()[i][d];
            if (this.cout(srce, dest) == 1 + this.cout(v, dest)) {
              this._tabRoutage[i][j].push(v);
            }
          }
        }
      }
    }
  };

  Eval.prototype.NCFT = function (srce, dest) {
    var som = 0;
    var v;
    if (srce.getId() == dest.getId()) return 0;
    else {
      if (this._graph.estVoisin(srce, dest) == 1) return 1;
      else {
        for (
          var d = 0;
          d < this._tabRoutage[srce.getId()][dest.getId()].length;
          d++
        ) {
          v = this._tabRoutage[srce.getId()][dest.getId()][d];
          som = som + this.NCFT(v, dest);
        }
        return som;
      }
    }
  };
  Eval.prototype.EstVoisinACO = function (entr, srce, dest) {
    for (var i in this._solACO) {
      if (
        this._solACO[i]._noeudEnt == entr.getId() &&
        this._solACO[i]._noeudSrc == srce.getId() &&
        this._solACO[i]._noeuddest == dest.getId()
      ) {
        return 1;
      }
    }
    return 0;
  };

  Eval.prototype.EstVoisinEuler = function (entr, srce, dest) {
    for (var i in this._solEuler) {
      if (
        this._solEuler[i]._noeudEnt == entr.getId() &&
        this._solEuler[i]._noeudSrc == srce.getId() &&
        this._solEuler[i]._noeuddest == dest.getId()
      ) {
        return 1;
      }
    }
    return 0;
  };

  Eval.prototype.NCCMACO = function (srce, dest) {
    var som = 0;
    var v;
    if (srce.getId() == dest.getId()) {
      return 0;
    } else {
      if (this._graph.estVoisin(srce, dest) == 1) return 1;
      else {
        for (
          var d = 0;
          d < this._tabRoutage[srce.getId()][dest.getId()].length;
          d++
        ) {
          v = this._tabRoutage[srce.getId()][dest.getId()][d];
          som = som + this.NCCMACO1(1, srce, v, dest);
        }
        return som;
      }
    }
  };
  Eval.prototype.NCCMEuler = function (srce, dest) {
    var som = 0;
    var v;
    if (srce.getId() == dest.getId()) {
      return 0;
    } else {
      if (this._graph.estVoisin(srce, dest) == 1) return 1;
      else {
        for (
          var d = 0;
          d < this._tabRoutage[srce.getId()][dest.getId()].length;
          d++
        ) {
          v = this._tabRoutage[srce.getId()][dest.getId()][d];
          som = som + this.NCCMEuler1(1, srce, v, dest);
        }
        return som;
      }
    }
  };

  Eval.prototype.NCCMACO1 = function (nivC, ne, s, d) {
    var som = 0;
    var v;
    if (s.getId() == d.getId()) {
      return 0;
    } else {
      if (
        this._graph.estVoisin(s, d) == 1 &&
        (this.EstVoisinACO(ne, s, d) == 1 || nivC < this._niveaux)
      )
        return 1;
      else {
        for (
          var m = 0;
          m < this._tabRoutage[s.getId()][d.getId()].length;
          m++
        ) {
          v = this._tabRoutage[s.getId()][d.getId()][m];
          if (this.EstVoisinACO(ne, s, v) == 1) {
            som = som + this.NCCMACO1(nivC, s, v, d);
          } else if (nivC < this._niveaux)
            som = som + this.NCCMACO1(nivC + 1, s, v, d);
        }
        return som;
      }
    }
  };

  Eval.prototype.NCCMEuler1 = function (nivC, ne, s, d) {
    var som = 0;
    var v;
    if (s.getId() == d.getId()) {
      return 0;
    } else {
      if (
        this._graph.estVoisin(s, d) == 1 &&
        (this.EstVoisinEuler(ne, s, d) == 1 || nivC < this._niveaux)
      )
        return 1;
      else {
        for (
          var m = 0;
          m < this._tabRoutage[s.getId()][d.getId()].length;
          m++
        ) {
          v = this._tabRoutage[s.getId()][d.getId()][m];
          if (this.EstVoisinEuler(ne, s, v) == 1) {
            som = som + this.NCCMEuler1(nivC, s, v, d);
          } else if (nivC < this._niveaux)
            som = som + this.NCCMEuler1(nivC + 1, s, v, d);
        }
        return som;
      }
    }
  };

  Eval.prototype.run = function (numEval) {
    var s, d;
    var ncft;
    var aco = 0.0;
    var euler = 0.0;
    for (var sind in this._graph.getNoeuds()) {
      s = this._graph.getNoeuds()[sind];
      for (var dind in this._graph.getNoeuds()) {
        d = this._graph.getNoeuds()[dind];
        if (s.getId() != d.getId()) {
          ncft = this.NCFT(s, d);
          if (numEval == 2) {
            aco = aco + this.NCCMACO(s, d) / ncft;
          }
          euler = euler + this.NCCMEuler(s, d) / ncft;
        }
      }
    }
    var nbNoeuds = this._graph.getNoeuds().length;
    if (numEval == 2) {
      aco = (aco * 100) / (nbNoeuds * (nbNoeuds - 1));
      this._metriqueACO = aco;
    }
    euler = (euler * 100) / (nbNoeuds * (nbNoeuds - 1));
    this._metriqueEuler = euler;
  };
  return Eval;
})();

var Tour = (function () {
  function Tour(graphA) {
    this._graphA = graphA;
    this._tour = [];
  }

  Tour.prototype.size = function () {
    return this._tour.length;
  };

  Tour.prototype.contains = function (noeudA) {
    for (var tourIndex in this._tour) {
      if (noeudA.isEqual(this._tour[tourIndex])) {
        return true;
      }
    }

    return false;
  };

  Tour.prototype.addNoeudA = function (noeudA) {
    this._tour.push(noeudA);
  };

  Tour.prototype.getNoeudA = function (tourIndex) {
    return this._tour[tourIndex];
  };

  return Tour;
})();

var AntColony = (function () {
  function AntColony(params, g, gA) {
    this._graph = g;
    this._graphA = gA;

    this._colony = [];

    this._colonySize = 5;
    this._alpha = 0.1;
    this._beta = 0.1;
    this._rho = 0.1;
    this._initPheromone = 1;
    this._maxIterations = 1;

    this.setParams(params);

    this._iteration = 0;

    this._tab = [];

    this._iterationBest = null;
    this._globalBest = null;

    this._createAnts();
  }

  AntColony.prototype.getGraph = function () {
    return this._graph;
  };
  AntColony.prototype.setGraph = function (g) {
    this._graph = g;
  };
  AntColony.prototype.getGraphA = function () {
    return this._graphA;
  };
  AntColony.prototype.setGraphA = function (gA) {
    this._graphA = gA;
  };

  AntColony.prototype.getColony = function () {
    return this._colony;
  };
  AntColony.prototype.size = function () {
    return this._colonySize;
  };
  AntColony.prototype.currentIteration = function () {
    return this._iteration;
  };
  AntColony.prototype.getMaxIterations = function () {
    return this._maxIterations;
  };
  AntColony.prototype.getiterationBest = function () {
    return this._iterationBest;
  };
  AntColony.prototype.getGlobalBest = function () {
    return this._globalBest;
  };
  AntColony.prototype.getTab = function () {
    return this._tab;
  };
  AntColony.prototype.getTabLigne = function (i) {
    return this._tab[i];
  };

  AntColony.prototype._createAnts = function () {
    this._colony = [];
    for (var antIndex = 0; antIndex < this._colonySize; antIndex++) {
      this._colony.push(
        new Ant(this._graph, this._graphA, {
          alpha: this._alpha,
          beta: this._beta,
          rho: this._rho,
        })
      );
    }
  };

  AntColony.prototype.setParams = function (params) {
    if (params != undefined) {
      if (params.colonySize != undefined) {
        this._colonySize = params.colonySize;
      }
      if (params.alpha != undefined) {
        this._alpha = params.alpha;
      }
      if (params.beta != undefined) {
        this._beta = params.beta;
      }
      if (params.rho != undefined) {
        this._rho = params.rho;
      }
      if (params.iteration != undefined) {
        this._maxIterations = params.iteration;
      }
      if (params.initPheromone != undefined) {
        this._initPheromone = params.initPheromone;
      }
    }
  };

  AntColony.prototype.reset = function () {
    this._iteration = 0;
    this._globalBest = null;
    this.resetAnts();
    this.setInitialPheromone(this._initPheromone);
    this._graph.resetPheromone();
  };

  AntColony.prototype.setInitialPheromone = function () {
    var edges = this._graph.getEdges();
    for (var edgeIndex in edges) {
      edges[edgeIndex].setInitialPheromone(this._initPheromone);
    }
  };

  AntColony.prototype.resetAnts = function () {
    this._createAnts();
    this._iterationBest = null;
  };

  AntColony.prototype.ready = function () {
    if (this._graph.size() <= 1) {
      return false;
    }
    return true;
  };

  AntColony.prototype.run = function () {
    if (!this.ready()) {
      return;
    }

    this._iteration = 0;
    while (this._iteration < this._maxIterations) {
      this.step();
    }
  };

  AntColony.prototype.algo = function () {
    var T_arc = this._graph.getTarc();
    var nbArc =
      2 *
      this._graph.getDimension() *
      Math.pow(2, this._graph.getDimension() - 1);
    var tabEuler;
    for (var t = 0; t < this._maxIterations; t++) {
      //pour chaque génération de fourmis
      tabEuler = this._graphA.cycleEuler();
      console.log("nb-arc= " + tabEuler.length);
      console.log("colony size = " + this._colonySize);
      this._createAnts();
      this._iterationBest = null;
      for (var n = 0; n < this._colonySize; n++) {
        //pour chaque fourmi
        console.log("fourmis num  " + n);
        //initialisation des colonnes 2 et 3 de T_arc
        for (var i = 0; i < nbArc; i++) {
          T_arc[i][2] = new Array();
          T_arc[i][3] = new Array();
        }
        this._colony[n].run(T_arc, tabEuler);
        //dépot de la phéromone
        this._colony[n].addPheromone();
        console.log(
          "la solution construise par la fourmis " +
            n +
            " a une longueur de " +
            this._colony[n].getTour().length
        );
      }
      this.getGlobalBest();
      //mise à jours de la phéromone
      for (var k = 0; k < this._graphA.getEdges().length; k++) {
        this._graphA
          .getEdges()
          [k].setPheromone(
            this._graphA.getEdges()[k].getPheromoneAvant() * this._rho +
              this._graphA.getEdges()[k].getPheromone() -
              this._graphA.getEdges()[k].getPheromoneAvant()
          );
        this._graphA
          .getEdges()
          [k].setPheromoneAvant(this._graphA.getEdges()[k].getPheromone());
      }
    }
  };

  AntColony.prototype.getIterationBest = function () {
    if (this._colony[0].getTour() == null) {
      return null;
    }

    if (this._iterationBest == null) {
      var best = this._colony[0];

      for (var antIndex in this._colony) {
        if (best.getTour().length < this._colony[antIndex].getTour().length) {
          best = this._colony[antIndex];
        }
      }
      this._iterationBest = best;
    }

    return this._iterationBest;
  };

  AntColony.prototype.getGlobalBest = function () {
    var bestAnt = this.getIterationBest();
    if (bestAnt == null && this._globalBest == null) {
      return null;
    }

    if (bestAnt != null) {
      if (
        this._globalBest == null ||
        this._globalBest.getTour().length <= bestAnt.getTour().length
      ) {
        this._globalBest = bestAnt;
      }
    }

    return this._globalBest;
  };

  return AntColony;
})();

var Ant = (function () {
  function Ant(graph, graphA, params) {
    this._graph = graph;
    this._graphA = graphA;
    this._alpha = params.alpha;
    this._beta = params.beta;
    this._rho = params.rho;
    this._tour = [];
  }

  Ant.prototype.reset = function () {
    this._tour = null;
  };

  Ant.prototype.getTour = function () {
    return this._tour;
  };

  Ant.prototype.calculProba = function (arc, tabEuler, sens) {
    var Tproba = [];
    var indice = 0;
    var edgetabeuler;
    var cumul = 0;
    var edge;
    var index;
    var ind;
    var heurist = 0.0;
    var max;
    var indMax;

    var tailleEuler = tabEuler.length;
    if (sens == 1) {
      index = this._graphA.getIndiceEuler(tabEuler, arc.getId());

      for (var k = 1; k <= tabEuler.length; k++) {
        index++;
        ind = index % tailleEuler;
        edge = this._graphA.getEdges()[tabEuler[ind]];
        if (edge.getVisite() == 0) {
          Tproba[indice] = new Array();
          Tproba[indice].push(edge.getId());
          Tproba[indice].push(1 / k);

          cumul =
            cumul +
            Math.pow(Tproba[indice][1], this._beta) *
              Math.pow(
                this._graphA.getEdges()[Tproba[indice][0]].getPheromone(),
                this._alpha
              );

          indice++;
        }
      }

      max = 0.0;
      for (var k = 0; k < Tproba.length; k++) {
        Tproba[k].push(
          (Math.pow(Tproba[k][1], this._beta) *
            Math.pow(
              this._graphA.getEdges()[Tproba[k][0]].getPheromone(),
              this._alpha
            )) /
            cumul
        );
        if (max < Tproba[k][2]) {
          max = Tproba[k][2];
          indMax = Tproba[k][0];
        }
      }
    } else if (sens == -1) {
      edgetabeuler = this._graphA.findEdge(arc.getDest(), arc.getSrc());

      index = this._graphA.getIndiceEuler(tabEuler, edgetabeuler.getId());

      for (var k = 1; k <= tabEuler.length; k++) {
        index--;
        ind = index % tailleEuler;

        if (index == -1) {
          index = tabEuler.length - 1;
        }
        edge = this._graphA.getEdges()[tabEuler[index]];
        edge = this._graphA.findEdge(edge.getDest(), edge.getSrc());

        if (edge.getVisite() == 0) {
          Tproba[indice] = new Array();
          Tproba[indice].push(edge.getId());
          Tproba[indice].push((1.0 / k) * tabEuler.length * 2);

          cumul =
            cumul +
            Math.pow(Tproba[indice][1], this._beta) *
              Math.pow(
                this._graphA.getEdges()[Tproba[indice][0]].getPheromone(),
                this._alpha
              );

          indice++;
        }
      }

      max = 0.0;
      var info = 0.0;
      var ph = 0.0;
      var res = 0.0;
      for (var k = 0; k < Tproba.length; k++) {
        info = Math.pow(Tproba[k][1], this._beta);
        ph = Math.pow(
          this._graphA.getEdges()[Tproba[k][0]].getPheromone(),
          this._alpha
        );
        res = (info * ph) / cumul;
        Tproba[k].push(res);
        if (max < Tproba[k][2]) {
          max = Tproba[k][2];
          indMax = Tproba[k][0];
        }
      }
    }

    return indMax;
  };

  Ant.prototype.run = function (T_arc, tabEuler) {
    var r;
    var arc;
    var index;
    var tour = [];
    //console.log("-------------sens direct------------");
    r = Math.floor(Math.random() * (tabEuler.length - 1));
    //console.log("l'indice choisi est "+r);
    arc = this._graphA.getEdges()[tabEuler[r]];
    this._graphA.getEdges()[arc.getId()].setVisite(1);
    //remplir les deux colonnes des pères et fils de arc dans T_arc
    T_arc = this._graphA.setPeresFils(tour, arc, T_arc);
    //parcours des arcs dans un seul sens
    for (var l = 0; l < this._graphA.getEdges().length / 2 - 1; l++) {
      // calcule de proba pour choisir le prochain saut
      arc = this._graphA.getEdges()[this.calculProba(arc, tabEuler, 1)];
      this._graphA.getEdges()[arc.getId()].setVisite(1);
      T_arc = this._graphA.setPeresFils(tour, arc, T_arc); //rajout à la solution inclu, si pas de cycle
    }
    //console.log("-------------sens inverse------------");
    arc = this._graphA.findEdge(arc.getDest(), arc.getSrc());
    T_arc = this._graphA.setPeresFils(tour, arc, T_arc);
    this._graphA.getEdges()[arc.getId()].setVisite(1);
    //parcours des arcs dans le sens inverse
    for (var l = 0; l < this._graphA.getEdges().length / 2 - 1; l++) {
      // calcule de proba pour choisir le prochain saut
      arc = this._graphA.getEdges()[this.calculProba(arc, tabEuler, -1)];
      this._graphA.getEdges()[arc.getId()].setVisite(1);
      T_arc = this._graphA.setPeresFils(tour, arc, T_arc); //rajout à la solution inclu, si pas de cycle
    }

    this._tour = tour;
    //réinitialiser le champs visité des arcs à 'non visité'
    this._graphA.resetVisite();
  };

  Ant.prototype.addPheromone = function () {
    var ph;
    for (var index in this._tour) {
      this._graphA
        .getEdges()
        [this._tour[index].getId()].addPheromone(Math.sqrt(this._tour.length));
    }
  };
  return Ant;
})();

var Graph = (function () {
  function Graph() {
    this._noeuds = [];
    this._tableVoisin = [];
    this._dimension = 6;
  }

  Graph.prototype.createG = function () {
    for (var i = 0; i < Math.pow(2, this._dimension); i++) {
      var noeudC = new Noeud(i);
      noeudC._createCoord(this._dimension);
      this._addNoeud(noeudC);
    }

    for (var NCindex = 0; NCindex < this.size(); NCindex++) {
      this._tableVoisin[NCindex] = new Array();

      console.log("////////////");

      for (var dim = 0; dim < this._dimension; dim++) {
        var d = new Noeud(-1);
        var table = Object.assign([], this.getNoeud(NCindex).getTableCoord());
        d.setTabCoor(table);

        d.setCoor(dim, 1 ^ this._noeuds[NCindex].getCoord(dim));
        d._calculId();

        this._tableVoisin[this._noeuds[NCindex].getId()][dim] = d;
      }
    }
  };

  Graph.prototype.estVoisin = function (s, d) {
    var compt = 0;
    for (var i = 0; i < this._dimension && compt < 2; i++) {
      if ((s.getCoord(i) ^ d.getCoord(i)) == 1) compt++;
    }
    if (compt == 1) return 1;
    else return 0;
  };

  Graph.prototype.getNoeuds = function () {
    return this._noeuds;
  };
  Graph.prototype.getTableVoisins = function () {
    return this._tableVoisin;
  };
  Graph.prototype.getVoisins = function (i) {
    return this._tableVoisin[i];
  };

  Graph.prototype.resetNoeudsVisit = function () {
    for (var ind in this._noeuds) {
      this._noeuds[ind].setVisit(0);
    }
  };

  // Graph.prototype.getBase = function() { return this._base; };
  Graph.prototype.getDimension = function () {
    return this._dimension;
  };

  Graph.prototype.getNoeud = function (noeudIndex) {
    return this._noeuds[noeudIndex];
  };

  Graph.prototype.size = function () {
    return this._noeuds.length;
  };

  Graph.prototype.setParam = function (params) {
    if (params != undefined) {
      //    if (params.base != undefined) {
      //       this._base = params.base;
      //    }
      if (params.dimension != undefined) {
        this._dimension = params.dimension;
      }
    }
  };

  Graph.prototype._addNoeud = function (noeud) {
    this._noeuds.push(noeud);
  };

  Graph.prototype.getTarc = function () {
    var T_arc = [];
    var nbArc = 2 * this._dimension * Math.pow(2, this._dimension);
    for (var i = 0; i < nbArc; i++) {
      T_arc[i] = new Array();
    }
    var count = 0;

    for (var i = 0; i < this._tableVoisin.length; i++) {
      for (var j = 0; j < this._dimension; j++) {
        T_arc[count].push(i);
        T_arc[count].push(this._tableVoisin[i][j].getId());
        count++;
      }
    }
    return T_arc;
  };

  Graph.prototype.resetPheromone = function () {
    for (var edgeIndex in this._edges) {
      this._edges[edgeIndex].resetPheromone();
    }
  };

  Graph.prototype.clear = function () {
    this._noeuds = [];
    this._tableVoisin = [];
    this._dimension = 1;
    this._base = 2;
  };

  return Graph;
})();

var Noeud = (function () {
  function Noeud(id) {
    this._id = id;
    this._tableCoord = [];
    this._visit = 0;
  }

  Noeud.prototype.getId = function () {
    return this._id;
  };
  Noeud.prototype.getTableCoord = function () {
    return this._tableCoord;
  };
  Noeud.prototype.getCoord = function (dimension) {
    return this._tableCoord[dimension];
  };
  Noeud.prototype.getCoorSize = function () {
    return this._tableCoord.length;
  };
  Noeud.prototype.getId = function () {
    return this._id;
  };
  Noeud.prototype.setVisit = function (val) {
    this._visit = val;
  };
  Noeud.prototype.setTabCoor = function (val) {
    this._tableCoord = val;
  };
  Noeud.prototype.setCoor = function (dim, coor) {
    this._tableCoord[dim] = coor;
  };

  Noeud.prototype._createCoord = function (dimension) {
    var id = this._id;
    this._tableCoord = new Array();
    for (var i = 0; i < dimension; i++) {
      this._addCoord(id % 2);
      id = Math.trunc(id / 2);
    }

    this._tableCoord.reverse();
  };

  Noeud.prototype._calculId = function () {
    var id = 0;
    var count = 0;
    for (var i = this._tableCoord.length - 1; i >= 0; i--) {
      id = id + this._tableCoord[count] * Math.pow(2, i);
      count++;
    }
    this._id = id;
  };

  Noeud.prototype._addCoord = function (coord) {
    this._tableCoord.push(coord);
  };

  return Noeud;
})();

Graph.prototype.createHgrapheA = function () {
  var gA = new GraphA();
  var TabDep = [];
  var nbarc = this._dimension * Math.pow(2, this._dimension - 1); // nombre d'arc de l'Hypercube(dimension) ( d * 2^(d-1) )
  for (var i = 0; i < nbarc; i++) {
    TabDep[i] = new Array();
  }
  var indexEdge = 0;
  // on rempli la table de dependances avec les noeuds de la table des voisins et leurs dependances
  for (var s = 0; s < this._tableVoisin.length; s++) {
    for (var d = 0; d < this._tableVoisin[s].length; d++) {
      if (s < this._tableVoisin[s][d].getId()) {
        var noeud = new NoeudA(
          s,
          this._tableVoisin[s][d].getId(),
          indexEdge,
          this._dimension
        );
        TabDep[indexEdge].push(noeud);
        indexEdge++;
      }
    }
  }

  // on aura jusqu'ici une table de dependances contenant nos d * 2^(d-1) nouveaux noeuds, en attendant d'ajouter celles de leurs Arcs
  console.log(indexEdge);
  console.log(TabDep);

  var myID = 0; // id des Arcs
  var l = 0;
  for (var ind = 0; ind < nbarc; ind++) {
    // pour chaque noeud dans la table de dependances on cherche les dependances qui restent avec les autres
    var l = 0; // emplacement de l'arc dans la table des edges d'un noeud
    for (var index = 0; index < nbarc; index++) {
      if (index != ind) {
        if (
          TabDep[ind][0].getSrc() == TabDep[index][0].getSrc() || // on test si ces derniers partage la meme source ou destination
          TabDep[ind][0].getSrc() == TabDep[index][0].getDest()
        ) {
          TabDep[ind].push(TabDep[index][0]);
          var dep = new Edge(
            myID,
            TabDep[ind][0].getId(),
            TabDep[index][0].getId(),
            TabDep[ind][0].getSrc()
          );
          gA.addEdge(dep); // on ajoute l'edge au graphe

          TabDep[ind][0]._addEdge(l, dep); // on ajoute l'edge au noeud
          TabDep[ind][0]._addEdgeEtiq(l, null);
          l++;
          myID++;
        }
        if (
          TabDep[ind][0].getDest() == TabDep[index][0].getSrc() || // on test si ces derniers partage la meme source ou destination
          TabDep[ind][0].getDest() == TabDep[index][0].getDest()
        ) {
          TabDep[ind].push(TabDep[index][0]);
          var dep = new Edge(
            myID,
            TabDep[ind][0].getId(),
            TabDep[index][0].getId(),
            TabDep[ind][0].getDest()
          );
          gA.addEdge(dep);
          TabDep[ind][0]._addEdge(l, dep);
          TabDep[ind][0]._addEdgeEtiq(l, null);
          l++;
          myID++;
        }
      }
    }
    gA.addNoeudA(TabDep[ind][0]); // ajout du noeud avec ses dependances au graphe adjoint
  }
  TabDep = null;
  //   setFormData({
  //     ...formData,
  //     nextStep: "block",
  //     loading: "none",
  //     valueEnd1: 40,
  //     nbrNA: gA._noeuds.length,
  //     nbrEA: gA._edges.length,
  //     hgraphA: gA,
  //     hgraph: this,
  //   });
  return gA;
};

var GraphA = (function () {
  function GraphA() {
    this._noeuds = [];
    this._edges = [];
  }

  GraphA.prototype.getEdges = function () {
    return this._edges;
  };
  GraphA.prototype.getNoeuds = function () {
    return this._noeuds;
  };

  GraphA.prototype.getNoeud = function (noeudIndex) {
    return this._noeuds[noeudIndex];
  };

  GraphA.prototype.size = function () {
    return this._noeuds.length;
  };

  GraphA.prototype.findNoeudA = function (id) {
    var found = 0;
    var i = 0;
    while (found != 1) {
      if (this._noeuds[i].getId() == id) {
        found = 1;
      }
      i++;
    }
    return this._noeuds[i - 1];
  };

  GraphA.prototype.getNoeudA = function (src, dest) {
    for (var index in this._noeuds) {
      if (
        this._noeuds[index].getSrc() == src &&
        this._noeuds[index].getDest()
      ) {
        return this._noeuds[index].getId();
      }
    }
    return -1;
  };

  GraphA.prototype.addNoeudA = function (noeud) {
    this._noeuds.push(noeud);
  };

  GraphA.prototype.addEdge = function (edge) {
    this._edges.push(edge);
  };

  GraphA.prototype.findEdge = function (src, dest) {
    for (var index in this._edges) {
      if (
        this._edges[index].getSrc() == src &&
        this._edges[index].getDest() == dest
      ) {
        return this._edges[index];
      }
    }
    return -1;
  };

  GraphA.prototype.resetPheromone = function () {
    for (var edgeIndex in this._edges) {
      this._edges[edgeIndex].resetPheromone();
    }
  };

  GraphA.prototype.calculEuler = function () {
    var r = Math.floor(Math.random() * (this._noeuds.length - 1));
    var racine = this._noeuds[r];
    racine.getTableEdges()[0][1] = 1;
    racine.incCPT();
    var noeudAvant = racine;
    var noeudCourant = this._noeuds[racine.getTableEdges()[0][0].getDest()];
    noeudCourant.setMsg(1);

    var etiq;
    var nbConst = 0;
    var link;
    var x;
    var y;
    var c;
    while (nbConst < this._noeuds.length) {
      if (noeudCourant.getMsg() == 1) {
        noeudCourant.incCPT();
        link = noeudCourant.find1(noeudCourant.getId(), noeudAvant.getId());
        noeudCourant.addEtiq(link, noeudCourant.getCPT());
        x = noeudCourant.canal_libre2(
          noeudCourant.getTableEdges()[link][0].getEtiq()
        );
        if (x != -1) {
          noeudCourant.incCPT();
          noeudCourant.addEtiq(x, noeudCourant.getCPT());
          noeudAvant = noeudCourant;
          noeudCourant = this._noeuds[
            noeudAvant.getTableEdges()[x][0].getDest()
          ];
          noeudCourant.setMsg(1);
        } else {
          noeudCourant.incCPV();
          y = noeudCourant.find2(noeudCourant.getCPV());
          noeudAvant = noeudCourant;
          noeudCourant = this._noeuds[
            noeudAvant.getTableEdges()[y][0].getDest()
          ];
          noeudCourant.setMsg(0);
        }
      } else {
        noeudCourant.incCPV();
        if (noeudCourant.getCPV() == noeudCourant.getTableEdges().length) {
          this._noeuds[noeudCourant.getId()] = noeudCourant;
          nbConst++;
        } else {
          x = noeudCourant.canal_libre2(
            noeudCourant.getTableEdges()[noeudCourant.find2(1)][0].getEtiq()
          );
          if (x != -1) {
            for (var j = noeudCourant.getCPT(); j >= 2; j--) {
              c = noeudCourant.find2(j);
              noeudCourant.addEtiq(
                c,
                noeudCourant.getTableEdges().length + j - noeudCourant.getCPT()
              );
            }
            noeudCourant.setCPT(2);
            noeudCourant.addEtiq(x, 2);
            noeudAvant = noeudCourant;
            noeudCourant = this._noeuds[
              noeudAvant.getTableEdges()[x][0].getDest()
            ];
            noeudCourant.setMsg(1);
          } else {
            noeudCourant.incCPV();
            link = noeudCourant.find1(noeudCourant.getId(), noeudAvant.getId());
            y = noeudCourant.find2(
              (noeudCourant.getTableEdges()[link][1] + 1) %
                (noeudCourant.getTableEdges().length + 1)
            );
            noeudAvant = noeudCourant;
            noeudCourant = this._noeuds[
              noeudAvant.getTableEdges()[y][0].getDest()
            ];
            noeudCourant.setMsg(0);
            if (noeudAvant.getCPV() == noeudAvant.getTableEdges().length) {
              this._noeuds[noeudAvant.getId()] = noeudAvant;
              nbConst++;
            }
          }
        }
      }
    }
  };

  GraphA.prototype.getIndiceEuler = function (tabEuler, val) {
    for (var index in tabEuler) {
      if (tabEuler[index] == val) {
        return index;
      }
    }
    return -1;
  };

  GraphA.prototype.cycleEuler = function () {
    var tabEuler = [];
    var r;
    var arc;
    var index;

    r = Math.floor(Math.random() * (this.getEdges().length - 1));
    arc = this.getEdges()[r];
    var noeudCourant = this._noeuds[arc.getDest()];

    arc = noeudCourant.getTableEdges()[noeudCourant.find2(1)][0];

    var num = noeudCourant.getNum(arc.getId());
    var node = this._noeuds[arc.getDest()];
    arc = node.getTableEdges()[node.find1(arc.getDest(), arc.getSrc())][0];
    tabEuler.push(arc.getId());
    var ind = 0;
    var k = 0;

    for (var i = 1; i < this._edges.length / 2; i++) {
      ind = noeudCourant.find2(
        (num + 1) % (noeudCourant.getTableEdges().length + 1)
      );
      arc = noeudCourant.getTableEdges()[ind][0];
      tabEuler.push(arc.getId());
      noeudCourant = this._noeuds[arc.getDest()];
      k = noeudCourant.find1(arc.getDest(), arc.getSrc());
      num = noeudCourant.getTableEdges()[k][1];
    }
    return tabEuler;
  };

  GraphA.prototype.newPush = function (tab, elem) {
    for (var index in tab) {
      if (tab[index] == elem) {
        return;
      }
    }
    tab.push(elem);
  };

  GraphA.prototype.rechLigneTarc = function (T_arc, src, dest) {
    for (var index in T_arc) {
      if (T_arc[index][0] == src && T_arc[index][1] == dest) {
        return index;
      }
    }
    return -1;
  };

  GraphA.prototype.estPere = function (T_arc, arc1, arc2) {
    var indiceArc2 = this.rechLigneTarc(T_arc, arc2[0], arc2[1]);
    var idArc1 = this.rechLigneTarc(T_arc, arc1[0], arc1[1]);
    var pereArc2 = T_arc[indiceArc2][3];
    for (var index in pereArc2) {
      if (pereArc2[index] == idArc1) {
        return 1;
      }
    }
    return 0;
  };

  GraphA.prototype.setPeresFils = function (tour, arc, T_arc) {
    var ae = [];
    var as = [];
    if (this._noeuds[arc.getSrc()].getSrc() == arc.getEtiq()) {
      ae.push(this._noeuds[arc.getSrc()].getDest());
      ae.push(this._noeuds[arc.getSrc()].getSrc());
    } else {
      ae.push(this._noeuds[arc.getSrc()].getSrc());
      ae.push(this._noeuds[arc.getSrc()].getDest());
    }

    if (this._noeuds[arc.getDest()].getDest() == arc.getEtiq()) {
      as.push(this._noeuds[arc.getDest()].getDest());
      as.push(this._noeuds[arc.getDest()].getSrc());
    } else {
      as.push(this._noeuds[arc.getDest()].getSrc());
      as.push(this._noeuds[arc.getDest()].getDest());
    }

    var indice1 = this.rechLigneTarc(T_arc, ae[0], ae[1]);

    var indice2 = this.rechLigneTarc(T_arc, as[0], as[1]);

    if (indice1 == undefined || indice2 == undefined) {
      alert("indice 1 = " + indice1 + "indice 2 = " + indice2);
    }

    if (this.estPere(T_arc, as, ae) == 0) {
      arc.setDep(ae[0], ae[1], as[1]);
      this.addSolution(tour, arc);

      this.newPush(T_arc[indice1][2], indice2);
      for (var k = 0; k < T_arc[indice2][2].length; k++) {
        this.newPush(T_arc[indice1][2], T_arc[indice2][2][k]);
      }

      this.newPush(T_arc[indice2][3], indice1);
      for (var k = 0; k < T_arc[indice1][3].length; k++) {
        this.newPush(T_arc[indice2][3], T_arc[indice1][3][k]);
      }
    }
    return T_arc;
  };

  GraphA.prototype.addSolution = function (tour, sol) {
    for (var index in tour) {
      if (tour[index].getId() == sol.getId()) {
        return;
      }
    }
    tour.push(sol);
  };

  GraphA.prototype.verifVesite = function () {
    for (var index in this.getEdges()) {
      if (this.getEdges()[index].getVisite() == 0) {
        return -1;
        console.log(
          "l'arc " + this.getEdges()[index].getId() + " n'a pas été visité"
        );
      }
    }
    return 1;
  };

  GraphA.prototype.verifSansCycle = function (tour, sens) {
    var ind = 1;
    if ((sens = 1)) {
      for (var index = 0; index < tour.length / 2; index++) {
        if (index == ind && tour[index].getEtiq() == tour[0].getEtiq()) {
          console.log(
            ".........interblocage index : " +
              index +
              ",etiq " +
              tour[index].getEtiq() +
              ", index+1 : 0, etiq " +
              tour[0].getEtiq()
          );
          return -1;
        } else if (
          tour[index].getEtiq() == tour[ind].getEtiq() &&
          index != ind
        ) {
          console.log(
            ".........interblocage index : " +
              index +
              ",etiq " +
              tour[index].getEtiq() +
              ", index+1 : " +
              ind +
              ", etiq " +
              tour[ind].getEtiq()
          );
          return -1;
        }
        if (ind < tour.length / 2 - 1) {
          ind++;
        }
      }
      return 1;
    } else if ((sens = -1)) {
      ind = tour.length / 2 + 1;
      for (var index = tour.length / 2; index < tour.length; index++) {
        if (
          index == ind &&
          tour[index].getEtiq() == tour[tour.length / 2].getEtiq()
        ) {
          console.log(
            ".........interblocage index : " +
              index +
              ",etiq " +
              tour[index].getEtiq() +
              ", index+1 : 0, etiq " +
              tour[tour.length / 2].getEtiq()
          );
          return -1;
        } else if (
          tour[index].getEtiq() == tour[ind].getEtiq() &&
          index != ind
        ) {
          console.log(
            ".........interblocage index : " +
              index +
              ",etiq " +
              tour[index].getEtiq() +
              ", index+1 : " +
              ind +
              ", etiq " +
              tour[ind].getEtiq()
          );
          return -1;
        }
        if (ind < tour.length - 1) {
          ind++;
        }
      }
      return 1;
    }
  };

  GraphA.prototype.resetVisite = function () {
    for (var index in this.getEdges()) {
      this.getEdges()[index].setVisite(0);
    }
  };

  GraphA.prototype.clear = function () {
    this._noeuds = [];
    this._edges = [];
  };

  return GraphA;
})();

var NoeudA = (function () {
  function NoeudA(s, d, i, dimension) {
    this._id = i;
    this._src = s;
    this._dest = d;
    this._tableEdges = [];
    for (var i = 0; i < 2 * (dimension - 1); i++) {
      this._tableEdges[i] = new Array();
    }

    this._cpt = 0;
    this._cpv = 0;
    this._satVerif = 0;
  }

  NoeudA.prototype.getSrc = function () {
    return this._src;
  };
  NoeudA.prototype.getDest = function () {
    return this._dest;
  };
  NoeudA.prototype.getId = function () {
    return this._id;
  };

  NoeudA.prototype.getTableEdges = function () {
    return this._tableEdges;
  };

  NoeudA.prototype.getEdge = function (index) {
    return this._tableEdges[index - 1];
  };

  NoeudA.prototype.incCPT = function () {
    return this._cpt++;
  };
  NoeudA.prototype.incCPV = function () {
    return this._cpv++;
  };
  NoeudA.prototype.setMsg = function (val) {
    return (this._satVerif = val);
  };
  NoeudA.prototype.getMsg = function () {
    return this._satVerif;
  };
  NoeudA.prototype.getCPT = function () {
    return this._cpt;
  };
  NoeudA.prototype.getCPV = function () {
    return this._cpv;
  };

  NoeudA.prototype.getNum = function (id) {
    var found = 0;
    var index = 0;
    while (found != 1 && index < this._tableEdges.length) {
      if (this._tableEdges[index][0].getId() == id) {
        found = 1;
        return this._tableEdges[index][1];
      }
      index++;
    }
    return -1;
  };

  NoeudA.prototype.setCPT = function (val) {
    this._cpt = val;
  };

  NoeudA.prototype._addEdge = function (indice, edge) {
    this._tableEdges[indice].push(edge);
  };
  NoeudA.prototype._addEdgeEtiq = function (indice, etiq) {
    this._tableEdges[indice].push(etiq);
  };

  NoeudA.prototype.addEtiq = function (indic, etiq) {
    this._tableEdges[indic][1] = etiq;
  };

  NoeudA.prototype.isEqual = function (noeudA) {
    if (this._id == noeudA.getId()) {
      return true;
    } else {
      return false;
    }
  };

  NoeudA.prototype.find1 = function (src, dest) {
    var found = 0;
    var index = 0;
    for (var index in this._tableEdges) {
      if (
        this._tableEdges[index][0].getSrc() == src &&
        this._tableEdges[index][0].getDest() == dest
      ) {
        return index;
      }
    }
    return -1;
  };

  NoeudA.prototype.find2 = function (val) {
    var found = 0;
    var index = 0;
    val = val % (this._tableEdges.length + 1);
    if (val == 0) {
      val = 1;
    }

    for (var index in this._tableEdges) {
      if (this._tableEdges[index][1] == val) {
        return index;
      }
    }
    return -1;
  };

  NoeudA.prototype.canal_libre1 = function () {
    var found = 0;
    var index = 0;
    while (found == 0 && index < this._tableEdges.length) {
      if (this._tableEdges[index][1] == null) {
        found = 1;
        return index;
      }
      index++;
    }
    return -1;
  };

  NoeudA.prototype.canal_libre2 = function (e) {
    var found = 0;
    var index = 0;
    while (found == 0 && index < this._tableEdges.length) {
      if (
        this._tableEdges[index][1] == null &&
        this._tableEdges[index][0].getEtiq() != e
      ) {
        found = 1;
        return index;
      }
      index++;
    }
    return -1;
  };

  return NoeudA;
})();

var Edge = (function () {
  function Edge(id, s, d, etiq) {
    this._id = id;
    this._src = s;
    this._dest = d;
    this._etiquette = etiq;
    this._pheromoneAvant = 1;
    this._pheromone = 1;
    this._visite = 0;
  }
  Edge.prototype.setDep = function (ne, s, d) {
    this._noeudEnt = ne;
    this._noeudSrc = s;
    this._noeuddest = d;
  };
  Edge.prototype.getId = function () {
    return this._id;
  };
  Edge.prototype.getSrc = function () {
    return this._src;
  };
  Edge.prototype.getDest = function () {
    return this._dest;
  };
  Edge.prototype.getEtiq = function () {
    return this._etiquette;
  };
  Edge.prototype.getPheromone = function () {
    return this._pheromone;
  };
  Edge.prototype.getVisite = function () {
    return this._visite;
  };

  Edge.prototype.setVisite = function (val) {
    this._visite = val;
  };

  Edge.prototype.contains = function (noeud) {
    if (this._src == noeud.getId()) {
      return true;
    }
    if (this._dest == noeud.getId()) {
      return true;
    }
    return false;
  };

  Edge.prototype.setInitialPheromone = function (pheromone) {
    this._initPheromone = pheromone;
  };

  Edge.prototype.setPheromone = function (pheromone) {
    this._pheromone = pheromone;
  };

  Edge.prototype.setPheromoneAvant = function (pheromoneA) {
    this._pheromoneAvant = pheromoneA;
  };

  Edge.prototype.getPheromoneAvant = function () {
    return this._pheromoneAvant;
  };

  Edge.prototype.addPheromone = function (pheromone) {
    this._pheromone = this._pheromone + pheromone;
  };
  Edge.prototype.resetPheromone = function () {
    this._pheromone = this._initPheromone;
  };

  return Edge;
})();

// INITIALISATION DES PRAMS

var hgraph = null;
var hgraphA = null;
var seuil = 85;

//

const execCreation = () => {
  hgraph = new Graph();
  hgraph.createG();
  console.log(hgraph);

  hgraphA = hgraph.createHgrapheA();

  console.log(hgraphA);

  hgraphA.calculEuler();
};

// ACO

// ant colony init

var colonySize = 10;
var alpha = 0.25;
var beta = 0.1;
var rho = 0.1;
var iteration = 10;
var initPheromone = 1;

var ant_colonyACO = null;

const execACO = () => {
  var params = {
    colonySize: colonySize,
    alpha: alpha,
    beta: beta,
    rho: rho,
    iteration: iteration,
    initPheromone: initPheromone,
  };

  var longuerState = 0;
  var eval1State = 0;

  var ant_colony = new AntColony(params, hgraph, hgraphA);
  console.log(ant_colony);
  ant_colony.algo();
  console.log(
    "longueur de la sol ACO =" + ant_colony._globalBest.getTour().length
  );

  longuerState = ant_colony._globalBest.getTour().length;
  console.log("");
  var heuristique = new EulerHeuristic(hgraph, hgraphA);
  heuristique.run();
  console.log("-----------Evaluation - niv 1---------");
  var niv = 1;
  var evaluation = new Eval(hgraph, heuristique.getTour(), niv);
  evaluation.addACOsol(ant_colony._globalBest.getTour());
  evaluation.calculTableRoutage();
  evaluation.run(2);

  console.log("métrique pour ACO niv 1 est =" + evaluation._metriqueACO + "%");
  eval1State = evaluation._metriqueACO.toFixed(2);

  //   var seuille = seuil;
  //   while (evaluation._metriqueACO < seuille) {
  //     console.log("Seuille non atteint");
  //     console.log("Passage au niveau suivant");
  //     niv++;
  //     console.log("-----------Evaluation - niv " + niv + "---------");
  //     var evaluation = new Eval(hgraph, heuristique.getTour(), niv);
  //     evaluation.addACOsol(ant_colony._globalBest.getTour());
  //     evaluation.calculTableRoutage();
  //     evaluation.run(2);
  //     console.log(
  //       "métrique pour ACO niv " + niv + " est =" + evaluation._metriqueACO + "%"
  //     );
  //   }
  //   console.log("Seuille atteint");
  //   console.log("");
  //   console.log(
  //     "Résultat : Seuille atteint avec ACO au niveau " +
  //       niv +
  //       " avec une métrique = " +
  //       evaluation._metriqueACO +
  //       "%" +
  //       "un seuil " +
  //       seuille
  //   );
  //   console.log("");
  ant_colonyACO = ant_colony;
  //   setFormData({
  //     ...formData,
  //     valueEnd1: evaluation._metriqueACO.toFixed(2),
  //     niveau: niv,
  //     aco_loading: "block",
  //     evaACO1: eval1State,
  //     longueurACO: longuerState,
  //     ant_colonyACO: ant_colony,
  //   });
};

// EULER

var heuristique = null;
var metrique_eul = 0;
const execEuler = () => {
  console.log("valuuuues");
  console.log(hgraph);
  console.log(hgraphA);

  heuristique = new EulerHeuristic(hgraph, hgraphA);
  heuristique.run();
  console.log("fin euler");
  console.log("Longueur de la solution Euler =" + heuristique.getTour().length);
  console.log("");
  console.log("-----------Evaluation - niv 1---------");
  var niv = 1;
  var evaluation = new Eval(hgraph, heuristique.getTour(), niv);
  evaluation.calculTableRoutage();
  evaluation.run(1);

  console.log(
    "métrique pour Euler niv 1 est =" + evaluation._metriqueEuler + "%"
  );
  var seuille = seuil;
  while (evaluation._metriqueEuler < seuille) {
    console.log("Seuille non atteint");
    console.log("Passage au niveau suivant");
    niv++;
    console.log("-----------Evaluation - niv " + niv + "---------");
    var evaluation = new Eval(hgraph, heuristique.getTour(), niv);
    evaluation.calculTableRoutage();
    evaluation.run(1);
    console.log(
      "métrique pour Euler niv " +
        niv +
        " est =" +
        evaluation._metriqueEuler +
        "%"
    );
  }
  console.log("Seuille atteint");
  console.log("");
  console.log(
    "Résultat : Seuille atteint avec Euler au niveau " +
      niv +
      " avec une métrique = " +
      evaluation._metriqueEuler +
      "%"
  );
  metrique_eul = evaluation._metriqueEuler;
  console.log("");
  //   setFormData({
  //     ...formData,
  //     valueEnd_euler: evaluation._metriqueEuler.toFixed(2),
  //     euler_loading: "block",
  //     heuristiqueE: heuristique,
  //     niveaueul: niv,
  //   });
};

// COMPARAISON EULER HEURSITIQUE DEULER

/// EVALUATION 1

var levelcomm = 2;

const execEval1 = () => {
  var mC = metrique_eul;
  var mcx = parseFloat(mC);
  console.log(mcx);

  if (heuristique == null) {
    alert("Veillez lancer l'heuristique d'euler d'abord ");
  } else {
    var niv = levelcomm;
    console.log("Patientez, évaluation en cours ...");

    console.log("");
    console.log("****Evaluation 1 : Heuristique d'Euler vs Cycle Eulérien***");
    var evaluation = new Eval(hgraph, heuristique.getTour(), niv);
    evaluation.calculTableRoutage();
    evaluation.run(1);
    console.log("fin eval");
    console.log("métrique cycle eulérien = " + mC + " %");
    console.log(
      "métrique heuristique d'euler = " + evaluation._metriqueEuler + " %"
    );
    if (mC > evaluation._metriqueEuler) {
      console.log(
        "Avec " +
          niv +
          " niveau(x) de communication, la méthode du cycle eulérien est meilleure que celle de l'heuristique d'euler pour cette instance. L'amélioration obtenue est de " +
          (mC - evaluation._metriqueEuler) +
          " %"
      );
      console.log(
        "Avec " +
          niv +
          " niveau(x) de communication, la méthode du cycle eulérien est meilleure que celle de l'heuristique d'euler pour cette instance. L'amélioration obtenue est de " +
          (mC - evaluation._metriqueEuler) +
          " %"
      );
    } else if (mC < evaluation._metriqueEuler) {
      console.log(
        "Avec " +
          niv +
          " niveau(x) de communication, la méthode de l'heuristique d'euler  est meilleure que celle du cycle eulérien pour cette instance. L'amélioration obtenue est de " +
          (evaluation._metriqueEuler - mC) +
          " %"
      );
      console.log(
        "Avec " +
          niv +
          " niveau(x) de communication, la méthode de l'heuristique d'euler  est meilleure que celle du cycle eulérien pour cette instance. L'amélioration obtenue est de " +
          (evaluation._metriqueEuler - mC) +
          " %"
      );
    } else {
      console.log(
        "Avec " +
          niv +
          " niveau(x) de communication,, les deux méthodes donnent la mème performance pour cette instance"
      );
      console.log(
        "Avec " +
          niv +
          " niveau(x) de communication,, les deux méthodes donnent la mème performance pour cette instance"
      );
    }

    console.log("");
    //   setFormData({
    //     ...formData,
    //     amelioration: Math.abs(evaluation._metriqueEuler - mC).toFixed(2),
    //     eval_loading: "block",
    //     aco_loading: "block",
    //     heur1: "Euler niv" + niv,
    //     heur2: "HeuristiqueEuler niv " + niv,
    //     valueEnd1: mcx.toFixed(2),
    //     valueEnd_euler: evaluation._metriqueEuler.toFixed(2),
    //     niveau: "",
    //     niveaueul: "",
    //   });
  }
};

// EVALUATION 2

var levelcomm2 = 2;

const execEval2 = () => {
  if (heuristique == null) {
    console.log("Veillez lancer l'heuristique d'euler d'abord ");
  } else if (ant_colonyACO == null) {
    console.log("Veillez lancer ACO d'abord ");
  } else {
    var niv = levelcomm2;
    console.log("Patientez, évaluation en cours ...");
    console.log("");
    console.log("****Evaluation 2 : ACO vs Heuristique d'Euler***");
    var evaluation = new Eval(hgraph, heuristique.getTour(), niv);
    evaluation.addACOsol(ant_colonyACO._globalBest.getTour());
    evaluation.calculTableRoutage();

    evaluation.run(2);
    console.log("fin eval");
    console.log("métrique pour ACO=" + evaluation._metriqueACO + "%");
    console.log("métrique pour Euler=" + evaluation._metriqueEuler + "%");
    if (evaluation._metriqueACO > evaluation._metriqueEuler) {
      console.log(
        "Avec " +
          niv +
          " niveau(x) de communication, la méthode ACO est meilleure que celle de l'heuristique d'euler pour cette instance. L'amélioration obtenue est de " +
          (evaluation._metriqueACO - evaluation._metriqueEuler) +
          " %"
      );
      console.log(
        "Avec " +
          niv +
          " niveau(x) de communication, la méthode ACO est meilleure que celle de l'heuristique d'euler pour cette instance. L'amélioration obtenue est de " +
          (evaluation._metriqueACO - evaluation._metriqueEuler) +
          " %"
      );
    } else if (evaluation._metriqueACO < evaluation._metriqueEuler) {
      console.log(
        "Avec " +
          niv +
          " niveau(x) de communication, la méthode de l'heuristique d'euler  est meilleure que ACO. L'amélioration obtenue est de " +
          (evaluation._metriqueEuler - evaluation._metriqueACO) +
          " %"
      );
      console.log(
        "Avec " +
          niv +
          " niveau(x) de communication, la méthode de l'heuristique d'euler  est meilleure que ACO. L'amélioration obtenue est de " +
          (evaluation._metriqueEuler - evaluation._metriqueACO) +
          " %"
      );
    } else {
      console.log(
        "Avec " +
          niv +
          " niveau(x) de communication,, les deux méthodes donnent la mème performance pour cette instance"
      );
      console.log(
        "Avec " +
          niv +
          " niveau(x) de communication,, les deux méthodes donnent la mème performance pour cette instance"
      );
    }
  }
  // setFormData({
  //   ...formData,
  //   amelioration: Math.abs(
  //     evaluation._metriqueEuler - evaluation._metriqueACO
  //   ).toFixed(2),
  //   eval_loading: "block",
  //   aco_loading: "block",
  //   heur1: "ACO niv " + niv,
  //   heur2: "HeuristiqueEuler niv " + niv,
  //   valueEnd1: evaluation._metriqueACO.toFixed(2),
  //   valueEnd_euler: evaluation._metriqueEuler.toFixed(2),
  //   niveau: "",
  //   niveaueul: "",
  // });
  console.log("");
};

execCreation();
execACO();
// execEuler();
// execEval1();
// execEval2();

// var heuristique = new EulerHeuristic(hgraph, hgraphA);
// heuristique.run();
// console.log("fin euler");

// var evaluation = new Eval(hgraph, heuristique.getTour(), 1);
// evaluation.calculTableRoutage();

// evaluation.run(1);
// console.log("fin eval1");
// console.log("métrique pour Euler=" + evaluation._metriqueEuler);

// var evaluation = new Eval(hgraph, heuristique.getTour(), 2);
// evaluation.calculTableRoutage();

// evaluation.run(1);
// console.log("fin eval2");
// console.log("métrique pour Euler=" + evaluation._metriqueEuler);

// var evaluation = new Eval(hgraph, heuristique.getTour(), 3);
// evaluation.calculTableRoutage();

// evaluation.run(1);
// console.log("fin eval3");
// console.log("métrique pour Euler=" + evaluation._metriqueEuler);

// var evaluation = new Eval(hgraph, heuristique.getTour(), 4);
// evaluation.calculTableRoutage();

// evaluation.run(1);
// console.log("fin eval4");
// console.log("métrique pour Euler=" + evaluation._metriqueEuler);

// var evaluation = new Eval(hgraph, heuristique.getTour(), 5);
// evaluation.calculTableRoutage();

// evaluation.run(1);
// console.log("fin eval5");
// console.log("métrique pour Euler=" + evaluation._metriqueEuler);

// var evaluation = new Eval(hgraph, heuristique.getTour(), 6);
// evaluation.calculTableRoutage();

// evaluation.run(1);
// console.log("fin eval6");
// console.log("métrique pour Euler=" + evaluation._metriqueEuler);
