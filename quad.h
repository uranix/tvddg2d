template<> struct quadrature<0> {
	static constexpr std::array<double, 1> s = {0.50000000000000000000, };
	static constexpr std::array<double, 1> w = {1.0000000000000000000, };
	static constexpr std::array<double, 1> mu = {1.0000000000000000000, };
	static constexpr std::array<std::array<double, 3>, 2> F2F = {
		1.0000000000000000000, 0, 0, 
		0, 0, 1.0000000000000000000, 
	};
};

constexpr std::array<double, 1> quadrature<0>::s;
constexpr std::array<double, 1> quadrature<0>::w;
constexpr std::array<double, 1> quadrature<0>::mu;
constexpr std::array<std::array<double, 3>, 2> quadrature<0>::F2F;

template<> struct quadrature<1> {
	static constexpr std::array<double, 2> s = {0.21132486540518711775, 0.78867513459481288225, };
	static constexpr std::array<double, 2> w = {0.5000000000000000000, 0.5000000000000000000, };
	static constexpr std::array<double, 2> mu = {-0.3660254037844386468, 1.366025403784438647, };
	static constexpr std::array<std::array<double, 4>, 3> F2F = {
		1.0000000000000000000, 0, 0, 0, 
		-0.3660254037844386468, 0.866025403784438647, 0.866025403784438647, -0.3660254037844386468, 
		0, 0, 0, 1.000000000000000000, 
	};
};

constexpr std::array<double, 2> quadrature<1>::s;
constexpr std::array<double, 2> quadrature<1>::w;
constexpr std::array<double, 2> quadrature<1>::mu;
constexpr std::array<std::array<double, 4>, 3> quadrature<1>::F2F;

template<> struct quadrature<2> {
	static constexpr std::array<double, 3> s = {0.11270166537925831148, 0.50000000000000000000, 0.88729833462074168852, };
	static constexpr std::array<double, 3> w = {0.277777777777777778, 0.444444444444444444, 0.277777777777777778, };
	static constexpr std::array<double, 3> mu = {0.187836108965430519, -0.666666666666666667, 1.478830557701236148, };
	static constexpr std::array<std::array<double, 5>, 4> F2F = {
		1.0000000000000000000, 0, 0, 0, 0, 
		-0.4788305577012361475, 1.07582870727983802, 0.573775310549246946, -0.358609569093279341, 0.187836108965430519, 
		0.187836108965430519, -0.35860956909327934, 0.57377531054924695, 1.07582870727983802, -0.47883055770123615, 
		0, 0, 0, 0, 1.0000000000000000, 
	};
};

constexpr std::array<double, 3> quadrature<2>::s;
constexpr std::array<double, 3> quadrature<2>::w;
constexpr std::array<double, 3> quadrature<2>::mu;
constexpr std::array<std::array<double, 5>, 4> quadrature<2>::F2F;

template<> struct quadrature<3> {
	static constexpr std::array<double, 4> s = {0.069431844202973712388, 0.33000947820757186760, 0.66999052179242813240, 0.93056815579702628761, };
	static constexpr std::array<double, 4> w = {0.17392742256872693, 0.32607257743127307, 0.32607257743127307, 0.173927422568726929, };
	static constexpr std::array<double, 4> mu = {-0.11391719628198993, 0.4007615203116504, -0.81363244948692726, 1.52678812545726679, };
	static constexpr std::array<std::array<double, 6>, 5> F2F = {
		1.0000000000000000000, 0, 0, 0, 0, 0, 
		-0.526788125457266787, 1.1590524262142826, 0.49403752802054745, -0.21435895436195619, 0.20197432186638285, -0.11391719628198993, 
		0.286844324029660474, -0.5315758353980039, 0.7447315113683434, 0.7447315113683434, -0.5315758353980039, 0.2868443240296605, 
		-0.113917196281989931, 0.2019743218663829, -0.2143589543619562, 0.4940375280205475, 1.1590524262142826, -0.5267881254572668, 
		0, 0, 0, 0, 0, 1.000000000000000, 
	};
};

constexpr std::array<double, 4> quadrature<3>::s;
constexpr std::array<double, 4> quadrature<3>::w;
constexpr std::array<double, 4> quadrature<3>::mu;
constexpr std::array<std::array<double, 6>, 5> quadrature<3>::F2F;

template<> struct quadrature<4> {
	static constexpr std::array<double, 5> s = {0.046910077030668003601, 0.23076534494715845448, 0.50000000000000000000, 0.76923465505284154552, 0.95308992296933199640, };
	static constexpr std::array<double, 5> w = {0.1184634425280945, 0.2393143352496832, 0.2844444444444444, 0.2393143352496832, 0.11846344252809454, };
	static constexpr std::array<double, 5> mu = {0.0763586617958129, -0.267941652223388, 0.533333333333333, -0.8931583920000717, 1.5514080490943130, };
	static constexpr std::array<std::array<double, 7>, 6> F2F = {
		1.0000000000000000000, 0, 0, 0, 0, 0, 0, 
		-0.551408049094313013, 1.200518144781687, 0.4596060639299722, -0.1713312197562774, 0.1169847996191190, -0.1307284012760011, 0.0763586617958129, 
		0.341750342905758725, -0.624281371610193, 0.822575656031950, 0.645245202549709, -0.327449694154935, 0.3337428547052848, -0.191582990427575, 
		-0.191582990427574609, 0.333742854705285, -0.327449694154935, 0.645245202549709, 0.822575656031950, -0.624281371610193, 0.341750342905759, 
		0.076358661795812900, -0.130728401276001, 0.116984799619119, -0.171331219756277, 0.459606063929972, 1.200518144781687, -0.551408049094313, 
		0, 0, 0, 0, 0, 0, 1.000000000000000, 
	};
};

constexpr std::array<double, 5> quadrature<4>::s;
constexpr std::array<double, 5> quadrature<4>::w;
constexpr std::array<double, 5> quadrature<4>::mu;
constexpr std::array<std::array<double, 7>, 6> quadrature<4>::F2F;

template<> struct quadrature<5> {
	static constexpr std::array<double, 6> s = {0.033765242898423986094, 0.16939530676686774317, 0.38069040695840154568, 0.61930959304159845432, 0.83060469323313225683, 0.96623475710157601391, };
	static constexpr std::array<double, 6> w = {0.0856622461895852, 0.180380786524069, 0.233956967286346, 0.2339569672863455, 0.180380786524069, 0.0856622461895852, };
	static constexpr std::array<double, 6> mu = {-0.054712724329266, 0.191800014038668, -0.379227702114614, 0.616930055430489, -0.940462843176349, 1.5656732001510719, };
	static constexpr std::array<std::array<double, 8>, 7> F2F = {
		1.0000000000000000000, 0, 0, 0, 0, 0, 0, 0, 
		-0.565673200151071933, 1.224169543733884, 0.441328821719275, -0.151797044205049, 0.0899371916948382, -0.075118590390608, 0.0918660019279970, -0.054712724329266, 
		0.374789643025276996, -0.67912284599003, 0.865170378724330, 0.599275195826888, -0.262791959278870, 0.197685743933936, -0.232093445950928, 0.137087289709402, 
		-0.242140412405211713, 0.41783647164643, -0.39336577668927, 0.717669717448056, 0.717669717448056, -0.39336577668927, 0.417836471646428, -0.24214041240521, 
		0.137087289709402042, -0.23209344595093, 0.19768574393394, -0.26279195927887, 0.599275195826888, 0.86517037872433, -0.67912284599003, 0.37478964302528, 
		-0.05471272432926591, 0.09186600192800, -0.07511859039061, 0.08993719169484, -0.15179704420505, 0.44132882171928, 1.22416954373388, -0.56567320015107, 
		0, 0, 0, 0, 0, 0, 0, 1.00000000000000, 
	};
};

constexpr std::array<double, 6> quadrature<5>::s;
constexpr std::array<double, 6> quadrature<5>::w;
constexpr std::array<double, 6> quadrature<5>::mu;
constexpr std::array<std::array<double, 8>, 7> quadrature<5>::F2F;
