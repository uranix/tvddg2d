#pragma once

#include "cell.h"

template<typename vars, int _p>
struct mesh {
	const double Lx, Ly;
	const int Nx, Ny;
	const double hx, hy;



	mesh(const double Lx, const double Ly, const int Nx, const int Ny)
		: Lx(Lx), Ly(Ly), Nx(Nx), Ny(Ny), hx(Lx/Nx), hy(Ly/Ny)
	{
		
	}
};
