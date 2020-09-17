import statistics
import numpy as np
import matplotlib.pyplot as plt
import random
import math
from sklearn.cluster import DBSCAN
import fun
import matplotlib.patches as patches
import Hybrid_R2

def levyheading():
    beta = 1.5
    num = math.gamma(1 + beta) * np.sin(np.pi * beta / 2)
    den = math.gamma((1 + beta) / 2) * beta * 2 ** ((beta - 1) / 2)
    sigma_u = (num / den) ** (1 / beta)
    k = 0  # random.randint(-50, 50)
    u = np.random.normal(k, sigma_u ** 2, 1)
    v = np.random.normal(k, 1, 1)
    step_size = abs(u[0] / (abs(v[0]) ** (1 / beta)))
    randnum = [-1, 1]
    theta = randnum[random.randint(0, 1)] * np.pi * random.random()
    # print(theta)
    if theta < 0:
        theta = 2 * np.pi + theta
    # print(step_size)
    if step_size > 20:
        print(step_size)
    return theta, step_size


def resheading(heading, repulsion, levy, ix, iy, iter):

    heading_factor = 0
    repulsion_factor = 0
    levy_factor = 1

    rescostheta = (heading_factor * math.cos(heading) + repulsion_factor * math.cos(repulsion) + levy_factor * math.cos(
        levy))
    ressintheta = (heading_factor * math.sin(heading) + repulsion_factor * math.sin(repulsion) + levy_factor * math.sin(
        levy))
    resheading = math.atan2(ressintheta, rescostheta)
    ix = ix + speed * rescostheta
    iy = iy + speed * ressintheta

    return ix, iy, resheading


def global_variable(distance, j, rpos, grp):
    costheta = 0
    sintheta = 0
    total_dist = 0
    for k in [h for h in grp if h != j]:
        # print(k)
        theta = math.atan2(rpos[j][1] - rpos[k][1], rpos[j][0] - rpos[k][0])
        costheta = (1 / distance[j][k]) * (math.cos(theta)) + costheta
        sintheta = (1 / distance[j][k]) * (math.sin(theta)) + sintheta
        total_dist = (1 / distance[j][k]) + total_dist

    costheta = ((costheta / total_dist))
    sintheta = ((sintheta / total_dist))
    restheta = math.atan2(sintheta, costheta)

    return restheta


def wallavoid(ix, iy, heading, nearobsx, nearobsy, dist, ox, oy):
    newheading = heading + np.pi / 2
    newix = ix + speed * math.cos(newheading)
    newiy = iy + speed * math.sin(newheading)
    a, b, c = fun.lineFromPoints(nearobsx, nearobsy)
    temp_dist = fun.perpendicular_distance(newix, newiy, a, b, c)
    newdist = np.zeros(len(oy)) + 10000

    if min(dist) > temp_dist:
        newheading = heading - np.pi / 2
        newix = ix + speed * math.cos(newheading)
        newiy = iy + speed * math.sin(newheading)

    a, b, c = fun.lineFromPoints(nearobsx, nearobsy)
    updated_temp_dist = fun.perpendicular_distance(newix, newiy, a, b, c)

    for xw in range(len(oy)):
        obsx = ox[xw]
        obsy = oy[xw]
        a, b, c = fun.lineFromPoints(obsx, obsy)
        newdist[xw] = fun.perpendicular_distance(newix, newiy, a, b, c)

    if updated_temp_dist > min(newdist):
        newheading = -heading
        newix = ix + speed * math.cos(newheading)
        newiy = iy + speed * math.sin(newheading)

    return newheading, newix, newiy

def levy_update(levy_previous, t, deltat):
    k = 0.1
    if deltat < t:
        levyupdate = speed * t - k * levy_previous
        if levyupdate < 0:
           __ , levyupdate = levyheading()
    else:
        levyupdate = speed * t + k * levy_previous

    return levyupdate

def wall_reflection(heading, wallnum, ix, iy):

    if wallnum == 1 or wallnum == 3:
        heading = -heading
    else:
        heading = np.pi - heading

    newix = ix + speed * math.cos(heading)
    newiy = iy + speed * math.sin(heading)

    return heading, newix, newiy

class motionplanning:

    def __init__(self, j, rx, ry):
        self.levy_heading = 0
        self.step_size = 0
        self.state = 3
        self.check = 0
        self.current_distance = 0
        self.ix, self.iy = rx[j], ry[j]
        self.heading, _ = levyheading()
        self.wallcheck = 0

    def update(self, ox, oy, rx, ry, j, nbot, iter, rpos, grp):

        ix = rx[j]
        iy = ry[j]

        dist = np.zeros(len(oy)) + 10000
        for xw in range(len(oy)):
            obsx = ox[xw]
            obsy = oy[xw]
            a, b, c = fun.lineFromPoints(obsx, obsy)
            dist[xw] = fun.perpendicular_distance(ix, iy, a, b, c)

            if min(dist) == dist[xw]:
                nearobsx = obsx
                nearobsy = obsy
                wallnum = xw

        self.levy_heading, self.step_size, self.current_distance = motionplanning.levystepsize(self)
        self.global_heading = motionplanning.globalheading(self, grp, rpos, j, iter)
        # if iter == 0:

        if min(dist) < 1 and self.check == 0:
            self.state = 1

        # if min(dist) < 1:
        #     self.state = 1
        #     self.check = 0
        # elif self.check > 10:
        #     self.state = 3
            # current_distance[j] = step_size[j]

        if self.state == 3:
            self.ix, self.iy, self.heading = resheading(self.heading, self.global_heading, self.levy_heading, ix, iy,
                                                        iter)
            self.current_distance = self.current_distance + speed

        # if self.wallcheck < 5 and self.wallcheck > 0:
        #     self.wallcheck = self.wallcheck + 1
        # else:
        #     self.wallcheck = 0
        if self.check > 0:
            self.check = self.check + 1
        if self.check > 5:
            self.check = 0

        if self.state == 1:
            if self.check == 0 and self.wallcheck == 0 or True:
                # self.heading, self.ix, self.iy = wallavoid(ix, iy, self.heading, nearobsx, nearobsy, dist, ox, oy)
                self.heading, self.ix, self.iy = wall_reflection(self.heading, wallnum, self.ix, self.iy)
                # levy_heading[j] = self.heading
                # self.current_distance = self.step_size
                self.levy_heading = self.heading
                self.current_distance = self.current_distance + speed
                self.state = 3
                self.wallcheck = 1
            # else:
            #     # self.heading = heading[j]
            #     self.ix = rx[j] + speed * math.cos(self.heading)
            #     self.iy = ry[j] + speed * math.sin(self.heading)
            #     self.current_distance = self.current_distance + speed
            # self.check = self.check + 1


    def levystepsize(self):
        if abs(self.current_distance - self.step_size) < 0.2:
            self.levy_heading, self.step_size = Hybrid_R2.levyheading()
            self.current_distance = 0

        return self.levy_heading, self.step_size, self.current_distance

    def globalheading(self, grp, rpos, j, iter):
        try:
            if iter != 0:
                for ki in range(len(grp)):
                    for kj in grp[ki][0]:
                        if kj == j:
                            globalheading = global_variable(distance, kj, rpos, grp[ki][0])

            return globalheading
        except UnboundLocalError or TypeError:
            return 0


def main():
    # calc potential field

    global distance
    global eff
    global bmap
    global speed
    # intialization
    show_animation = True
    nbot = 10 # number of robots
    area = 40
    pixel = 100
    bmap = np.zeros([area*pixel, area*pixel])
    speed = 0.178

    oy = [[0, 0], [0, area], [area, area], [area, 0]]
    ox = [[0, area], [area, area], [area, 0], [0, 0]]
    rx = []
    ry = []

    for i in range(int(nbot)):
    #     rx.append((area * random.random()))
    #     ry.append((area * random.random()))
        # rx.append((3))
        # ry.append((area/2))
        rx.append((10 + 2 *random.random()))
        ry.append((2 + 1 *random.random()))
    # rx = [18.484494338927238, 24.4720554824703, 17.060682131639552, 30.02291051004108, 3.0048652509077582,
    #       37.16401227845603, 11.480780226537842, 21.12814300744732, 6.065440706295824, 34.73113424594841]
    # ry = [4.091447935267879, 20.317942519639217, 30.305645541495892, 38.40802245438094, 33.48078494040231,
    #       12.297637558915177, 12.240438073091664, 32.80082006697463, 12.956461388511034, 8.64519276887101]


    deltat = np.zeros(nbot)
    bot = np.zeros(nbot)
    bot = list(bot)
    botnew = np.zeros(nbot)
    botnew = list(botnew)
    iterations = 1500
    eff = np.zeros(iterations)
    nearbot = np.zeros([nbot]) + 100
    flagtemp = np.zeros([nbot])
    rpos = np.zeros([nbot, 2])
    distance = np.zeros([nbot, nbot])
    if show_animation:
        # rect = patches.Rectangle((1, 1), 18, 18, linewidth=0.5, edgecolor='r', facecolor='none')

        # Add the patch to the Axes
        # fig, ax = plt.subplots(1)
        plt.ylim(-1, 41)
        plt.xlim(-1, 41)
        # ax.imshow()
        # ax.add_patch(rect)

    for ri in range(nbot):
        for rj in range(nbot):
            distance[ri][rj] = np.hypot(ry[ri] - ry[rj], rx[ri] - rx[rj])
    grp = 0
    for j in range(nbot):
        bot[j] = motionplanning(j, rx, ry)
    for j in range(nbot):
        botnew[j] = motionplanning(j, rx, ry)
    for itn in range(iterations):
        # print(i)
        for j in range(nbot):
            bot[j].update(ox, oy, rx, ry, j, nbot, itn, rpos, grp)
            print(np.hypot(ry[j] - bot[j].iy, rx[j] - bot[j].ix))
            rx[j] = bot[j].ix
            ry[j] = bot[j].iy

            # updating efficiency

            z = (np.linspace(-30, 30, 61))
            for ix in z:
                for iy in z:
                    try:
                        bmap[int(rx[j]*pixel) + int(ix), int(ry[j]*pixel) + int(iy)] = 1
                    except IndexError:
                        pass
            # rmap[j][int(rx[j])][int(ry[j])] = 1

        eff[itn] = (np.count_nonzero(bmap == 1))

        for w in range(nbot):
            rpos[w][0] = bot[w].ix
            rpos[w][1] = bot[w].iy

        # clustering for DBSCAN
        clustering = DBSCAN(eps=2, min_samples=2).fit(rpos)
        grp = []
        numclus = max(clustering.labels_)
        for w1 in range(numclus + 1):
            grp.append(np.where(clustering.labels_ == w1))


        # updating distance and delta t

        for ri in range(nbot):
            for rj in range(nbot):
                distance[ri][rj] = np.hypot(bot[ri].ix - bot[rj].ix, bot[ri].iy - bot[rj].iy)
        colors = [".r"]
        for ri in range(nbot):
            distance[ri][ri]=1000
            mindistance = min(distance[ri])
            rj = np.where(distance[ri] == mindistance)
            rj = rj[0][0]
            if mindistance < 0.6:
                meant = statistics.mean(deltat)
                if itn < 10:
                    colors = [".b"]
                    botnew[ri].step_size = 5
                    botnew[ri].levy_heading = global_variable(distance, ri, rpos, grp[0][0])
                else:
                    botnew[ri].step_size = levy_update(bot[ri].step_size, meant, deltat[ri])
                    botnew[ri].levy_heading = (math.atan2(bot[ri].iy - bot[rj].iy, bot[ri].ix - bot[rj].ix))

                deltat[ri] = 0
                nearbot[ri] = rj
                flagtemp[ri] = 0
            else:
                botnew[ri].step_size = bot[ri].step_size
                botnew[ri].levy_heading = bot[ri].levy_heading
            if flagtemp[ri] > 10:
                nearbot[ri] = 100

        flagtemp = flagtemp + 1
        deltat = deltat + 1

        for i in range(nbot):
            bot[i].step_size = botnew[i].step_size
            bot[i].levy_heading = botnew[i].levy_heading

        # plotting
        if show_animation:
            for j in range(nbot):
                # plt.ylim(0, 20)
                # plt.xlim(0, 20)

                # if bot[j].state == 1:
                    # plt.plot(bot[j].ix, bot[j].iy, colors[0])
                # if bot[j].state == 3:
                plt.plot(bot[j].ix, bot[j].iy, '.')
                # if j == 0:
                #     plt.pause(0.005)

    print("Goal!!")

    if show_animation:
        plt.show()

    plt.figure()
    draw_heatmap(bmap)
    return rx, ry


def draw_heatmap(data):
    data = np.array(data).T
    z_min, z_max = data.min(), data.max()
    plt.pcolormesh(data, vmin=z_min, vmax=z_max, cmap='RdBu')
    plt.show



if __name__ == '__main__':
    iterations = 600
    times = 1
    effmat = np.zeros((iterations, times))
    rate_inc = np.zeros((iterations, times))
    for i in range(times):
        print(__file__ + " start!!")
        main()
        bmap
        distance
        print(__file__ + " Done!!")
        effmat[:, i] = eff
        # plt.figure()
        # plt.plot(eff)
        # plt.show()

    # for k in range(1, iterations):
    #     rate_inc[k - 1] = (effmat[k] - effmat[k - 1]) * 100 / effmat[k - 1]
    #
    # plt.figure()
    # plt.plot(rate_inc)
    # plt.show()
