#include <cstdlib>
#include <iostream>

long hilbidx(long index, long sidePow) {
    long x = 0;
    long y = 0;
    long t = index;
    long s = 1;

    for (long i = 0; i < sidePow; i++) {
        long rx = 1 & (t / 2);
        long ry = 1 & (t ^ rx);

        if (ry == 0) {
            if (rx == 1) {
                x = s - 1 - x;
                y =  s - 1 - y;
            }

            long temp = x;
            x = y;
            y = temp;
        }

        x += s * rx;
        y += s * ry;
        t /= 4;
        s *= 2;
    }

    return x + (y << sidePow);
}

static int sign(int x) {
  if (x < 0) { return -1; }
  if (x > 0) { return  1; }
  return 0;
}

int in_bounds(int x,  int y,
               int x_s,int y_s,
               int ax, int ay,
               int bx, int by) {
    int dx, dy;

    dx = ax + bx;
    dy = ay + by;

    if (dx < 0) {
        if ((x > x_s) || (x <= (x_s + dx))) { return 0; }
    }
    else {
        if ((x < x_s) || (x >= (x_s + dx))) { return 0; }
    }

    if (dy < 0) {
        if ((y > y_s) || (y <= (y_s + dy))) { return 0; }
    }
    else {
        if ((y < y_s) || (y >= (y_s + dy))) { return 0; }
    }

    return 1;
}

int recgilbert(int cur_idx,
                   int x_dst, int y_dst,
                   int x, int y,
                   int ax, int ay,
                   int bx,int by ) {

    int width = std::abs(ax + ay);
    int height = std::abs(bx + by);
  
    // unit major direction
    int dax = sign(ax);
    int day = sign(ay);

    // unit orthogonal direction
    int dbx = sign(bx);
    int dby = sign(by);

    int dx = dax + dbx;
    int dy = day + dby;

    if (height == 1) {
        if (dax == 0) return cur_idx + dy * (y_dst - y);
        return cur_idx + dx * (x_dst - x);
    }

    if (width == 1) {
        if (dbx == 0) return cur_idx + dy * (y_dst - y);
        return cur_idx + dx * (x_dst - x);
    }

    int ax2 = ax >> 1;
    int ay2 = ay >> 1;
    int bx2 = bx >> 1;
    int by2 = by >> 1;

    int w2 = abs(ax2 + ay2);
    int h2 = abs(bx2 + by2);

    if (2 * width > 3 * height) {
        if ((w2 & 1) && (width > 2)) {
            // prefer even steps
            ax2 += dax;
            ay2 += day;
        }
        if (in_bounds(x_dst, y_dst, x, y, ax2, ay2, bx, by)) {
            return recgilbert(cur_idx, x_dst, y_dst, x, y, ax2, ay2, bx, by);
        }

        cur_idx += abs((ax2 + ay2) * (bx + by));
        return recgilbert(cur_idx, x_dst, y_dst, x+ax2, y+ay2, ax-ax2, ay-ay2, bx, by);
    }

    if ((h2 & 1) && (height > 2)) {
        // prefer even steps
        bx2 += dbx;
        by2 += dby;
    }

    // standard case, one step up, one long horizontal, one step down
    if (in_bounds(x_dst, y_dst, x, y, bx2, by2, ax2, ay2)) {
        return recgilbert(cur_idx, x_dst, y_dst, x, y, bx2, by2, ax2, ay2);
    }
    cur_idx += abs((bx2 + by2) * (ax2 + ay2));

    if (in_bounds(x_dst, y_dst, x+bx2, y+by2, ax, ay, bx-bx2, by-by2)) {
        return recgilbert(cur_idx, x_dst, y_dst, x+bx2, y+by2, ax, ay, bx-bx2, by-by2);
    }
    cur_idx += abs((ax + ay) * (bx - bx2 + by - by2));

    return recgilbert(cur_idx, x_dst, y_dst, x+ax-dax+bx2-dbx, y+ay-day+by2-dby, -bx2, -by2, -(ax-ax2), -(ay-ay2));
}

long gilbidx(long x, long y, long width, long height) {
    if (width >= height) return recgilbert(0, x, y, 0, 0, width, 0, 0, height);
    else return recgilbert(0, x, y, 0, 0, 0, height, width, 0);
}

// https://github.com/jakubcerveny/gilbert/blob/master/ports/gilbert.c