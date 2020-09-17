import statistics
import numpy as np
# import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
from celluloid import Camera
import random
import math
from sklearn.cluster import DBSCAN
import fun
import functions_movement
import matplotlib.patches as patches


def resheading(heading, repulsion, levy, ix, iy, iter):

    heading_factor = 0.6
    repulsion_factor = 0
    levy_factor = 0
    if repulsion == 0:
        repulsion_factor = 0

    rescostheta = (heading_factor * math.cos(heading) + repulsion_factor * math.cos(repulsion) + levy_factor * math.cos(
        levy))
    ressintheta = (heading_factor * math.sin(heading) + repulsion_factor * math.sin(repulsion) + levy_factor * math.sin(
        levy))
    resheading = math.atan2(ressintheta, rescostheta)
    ix = ix + speed * math.cos(resheading)
    iy = iy + speed * math.sin(resheading)
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
        self.heading = 0
        self.wallcheck = 0
        self.left_sensor, self.front_sensor, self.right_sensor = 0, 0, 0
        self.sensor_range = 0.6
        self.velocity = [0, 0]
        self.lbest_velocity = [0, 0]
        self.gbest_velocity = [0, 0]
        self.position = [rx[j], ry[j]]
        self.start = [rx[j], ry[j]]
        self.levy_1 = [rx[j], ry[j]]
        self.levy_2 = [rx[j], ry[j]]
        self.sumspeed = 0
        self.interactions = []

    def update(self, ox, oy, rx, ry, j, nbot, iter, rpos, grp, current_grid, pixel):

        ix = rx[j]
        iy = ry[j]

        self.position = [rx[j], ry[j]]

        self.left_sensor, self.front_sensor, self.right_sensor = functions_movement.sensor_reading(current_grid, self,
                                                                                                   pixel,
                                                                                                   self.sensor_range)

        if self.left_sensor + self.front_sensor + self.right_sensor > 0:
            self.state = 1
        else:
            self.state = 3

        if self.state == 1:
            self.velocity = functions_movement.obstacle_avoidance_velocity(self.left_sensor, self.front_sensor, self.right_sensor, self)
            self.step_size = 0
        elif iter%5 == 0:
            if (-self.current_distance + self.step_size) < 0.05:
                self.lbest_velocity, self.step_size = functions_movement.levy_movement()
                self.levy_2 = self.levy_1
                self.levy_1 = self.position
                self.current_distance = 0
            if random.random() < 0.02:
                self.lbest_velocity = motionplanning.pointrepel(self)
            #     print(j)


            self.gbest_velocity = motionplanning.globalheading(self, grp, rpos, j, iter)
            if np.hypot(self.gbest_velocity[1],self.gbest_velocity[0]) > 0:
            # print(self.gbest_velocity)
                inertia_factor = 0.6
                lbest_factor = 0.1 * random.random()
                gbest_factor = 2 * random.random()
            else:
                inertia_factor = 0.6
                lbest_factor = 2 * random.random()
                gbest_factor = 2 * random.random()

            for i in range(2):
                self.velocity[i] = inertia_factor * self.velocity[i] + lbest_factor * self.lbest_velocity[i]  + gbest_factor * self.gbest_velocity[i]

        self.heading = math.atan2(self.velocity[1], self.velocity[0])

        self.velocity[0] = self.velocity[0]
        self.velocity[1] = self.velocity[1]

        if np.hypot(self.velocity[1],self.velocity[0]) > 0.2:
            self.velocity[0] = 0.2*math.cos(self.heading)
            self.velocity[1] = 0.2*math.sin(self.heading)

        for i in range(2):
            self.position[i] = self.position[i] + self.velocity[i]

        self.sumspeed = self.sumspeed + np.hypot(self.velocity[1], self.velocity[0])

        self.ix = self.position[0]
        self.iy = self.position[1]
        self.heading = math.atan2(self.velocity[1], self.velocity[0])
        self.current_distance = self.current_distance + np.hypot(self.velocity[0], self.velocity[1])
        # print(np.hypot(self.velocity[0], self.velocity[1]))


    def globalheading(self, grp, rpos, j, iter):
        try:
            if iter != 0:
                for ki in range(len(grp)):
                    for kj in grp[ki][0]:
                        if kj == j:
                            gbest_velocity = functions_movement.global_variable(distance, kj, rpos, grp[ki][0])

            return gbest_velocity
        except UnboundLocalError or TypeError:
            return [0, 0]

    def pointrepel(self):
        velocity = [0, 0]
        if random.random() < 0.99:
            repel_point = self.levy_2
            print("im here idiot")
        else:
            repel_point = self.start

        start_distance = np.hypot(self.position[0]-repel_point[0], self.position[1]-repel_point[1])
        theta = math.atan2(self.position[1] - repel_point[1], self.position[0] - repel_point[0])

        if start_distance <0.00001:
            print("stop")
            start_distance = 0.1
        total_dist = (1 / start_distance)

        velocity[0] = total_dist * math.cos(theta)
        velocity[1] = total_dist * math.sin(theta)

        return velocity



def main():
    # calc potential field

    global distance
    global eff
    global bmap
    global speed
    global victim_found
    global bmap
    global camera
    # intialization
    show_animation = True
    nbot = 10  # number of robots
    nvictim = 4
    area = 40
    pixel = 100
    bmap = np.zeros([area*pixel, area*pixel])
    speed = 0.1
    fig = plt.figure()
    camera = Camera(fig)

    oy = [[0, 0], [0, area], [area, area], [area, 0]]
    ox = [[0, area], [area, area], [area, 0], [0, 0]]
    rx = []
    ry = []

    victims_x = []
    victims_y = []
    victim_found = np.zeros([nvictim])
    avg_speed = np.zeros([nbot])
    for i in range(int(nbot)):
        rx.append(area * random.random())
        ry.append(area * random.random())
        # rx.append((10))
        # ry.append(10)
        # rx.append((20 + 0.5 *random.random()))
        # ry.append((20 + 0.5 *random.random()))
    # rx = [18.484494338927238, 24.4720554824703, 17.060682131639552, 30.02291051004108, 3.0048652509077582, 37.16401227845603, 11.480780226537842, 21.12814300744732, 6.065440706295824, 34.73113424594841]
    # ry = [4.091447935267879,20.317942519639217,30.305645541495892,38.40802245438094,33.48078494040231,12.297637558915177,12.240438073091664,32.80082006697463,12.956461388511034, 8.64519276887101]
    for i in range(nvictim):
        victims_x.append((area * random.random()))
        victims_y.append((area * random.random()))

    plt.plot(victims_x, victims_y, '.r')

    bot = np.zeros(nbot)
    bot = list(bot)
    iter = 1500
    eff = np.zeros(iter)
    rpos = np.zeros([nbot, 2])
    distance = np.zeros([nbot, nbot])
    vic_distance = np.zeros([nvictim, nbot])
    if show_animation:
        # rect = patches.Rectangle((1, 1), 18, 18, linewidth=0.5, edgecolor='r', facecolor='none')

        # Add the patch to the Axes
        # fig, ax = plt.subplots(1)
        plt.ylim(-1, 41)
        plt.xlim(-1, 41)
        # ax.imshow()
        # ax.add_patch(rect)

    for ri in range(nbot):
        for rj in [h for h in range(nbot) if ri != h]:
            distance[ri][rj] = np.hypot(ry[ri] - ry[rj], rx[ri] - rx[rj])
    grp = 0
    for j in range(nbot):
        bot[j] = motionplanning(j, rx, ry)
    for itn in range(iter):
        # print("iteration num:")
        # print(itn)
        current_grid = np.zeros([area*pixel, area*pixel])
        current_grid[0:1,:]=1
        current_grid[area*pixel-2:area*pixel-1, :] = 1
        current_grid[:, 0:1] = 1
        current_grid[:, area * pixel - 2:area * pixel-1] = 1

        # for j in range(nbot):
        #     current_grid[int(rx[j]*pixel), int(ry[j]*pixel)] = 1

        for j in range(nbot):
            bot[j].update(ox, oy, rx, ry, j, nbot, itn, rpos, grp,current_grid, pixel)
            # print(np.hypot(ry[j] - bot[j].iy, rx[j] - bot[j].ix))
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
            for rj in [h for h in range(nbot) if ri != h]:
                distance[ri][rj] = np.hypot(bot[ri].ix - bot[rj].ix, bot[ri].iy - bot[rj].iy)
                if distance[ri][rj] < 2:
                    temp_interlen = len(bot[ri].interactions)
                    if temp_interlen > 0:
                        if bot[ri].interactions[temp_interlen-1][0] == rj:
                                if (itn - bot[ri].interactions[temp_interlen-1][1]) > 10:
                                    bot[ri].interactions.append([rj, itn])
                        else:
                            bot[ri].interactions.append([rj, itn])
                    else:
                        bot[ri].interactions.append([rj, itn])


        for ri in range(nvictim):
            for rj in range(nbot):
                vic_distance[ri][rj] = np.hypot(victims_x[ri] - bot[rj].ix, victims_y[ri] - bot[rj].iy)
        # plotting
        for ri in range(nvictim):
            if min(vic_distance[ri]) < 1:
                victim_found[ri] = 1
                # print(victim_found)
        # if sum(victim_found) == nvictim:
        #     break
        if show_animation:
            for j in range(nbot):

                # if bot[j].state == 1:
                plt.plot(bot[j].ix, bot[j].iy, '.b')

                # if bot[j].state == 3:
                #     plt.plot(bot[j].ix, bot[j].iy, '.b')
                # if j == 0:
                #     plt.pause(0.01)
        # camera.snap()

        for w1 in range(len(grp)):
            if len(grp[w1][0])>2:
                for j in range(nbot):
                    plt.plot(bot[j].ix, bot[j].iy, '.r', markersize=30)
                # plt.pause(0.01)
                print("cluster")



    # animation.save('celluloid_minimal.gif', writer='imagemagick')
    for i in range(nbot):
        avg_speed[i] = bot[i].sumspeed/(itn+1)

# def draw_heatmap(data):
#     data = np.array(data).T
#     z_min, z_max = data.min(), data.max()
#     plt.pcolormesh(data, vmin=z_min, vmax=z_max, cmap='RdBu')
#     plt.show



if __name__ == '__main__':
    global beta
    # global camera
    iterations = 600
    times = 1
    effmat = np.zeros((iterations, times))
    rate_inc = np.zeros((iterations, times))
    betamat = np.linspace(1.5, 1.5, num = 1)
    beta = 1.5
    for i in range(times):

        print(__file__ + " start!!")
        # beta = betamat[i]
        main()
        victim_found
        # bmap
        # distance
        print(__file__ + " Done!!")

        # effmat[:, i] = eff
    # draw_heatmap(bmap)
    # for k in range(1, iterations):
    #     rate_inc[k - 1] = (effmat[k] - effmat[k - 1]) * 100 / effmat[k - 1]

    # plt.figure()
    # plt.plot(rate_inc)
    # plt.show()
    # camera.animate()