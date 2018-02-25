#pragma once

struct grid {
    const int Nx, Ny;
    const double Lx, Ly;
    const double hx, hy;
    grid(const int Nx, const int Ny, const double Lx, const double Ly)
        : Nx(Nx), Ny(Ny), Lx(Lx), Ly(Ly), hx(Lx/Nx), hy(Ly/Ny)
    { }
};
