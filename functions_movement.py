import numpy as np
import math
import random


def sensor_reading(current_grid, bot, pixel, sensor_range):
    front_sensor = 0
    right_sensor = 0
    left_sensor = 0

    x = int(bot.ix * pixel)
    y = int(bot.iy * pixel)
    # current_grid[x, y] = 0
    theta = bot.heading
    for i in range(int(sensor_range * pixel)):
        try:
            if current_grid[int(x + i * math.cos(theta)), int(y + i * math.sin(theta))] > 0:
                front_sensor = 1
        except IndexError:
            pass
        try:
            if current_grid[int(x + i * math.cos(theta + np.pi / 2)), int(y + i * math.sin(theta + np.pi / 2))] > 0:
                left_sensor = 1
        except IndexError:
            pass
        try:
            if current_grid[int(x + i * math.cos(theta - np.pi / 2)), int(y + i * math.sin(theta - np.pi / 2))] > 0:
                right_sensor = 1
        except IndexError:
            pass

    return left_sensor, front_sensor, right_sensor


def obstacle_avoidance(left_sensor, front_sensor, right_sensor, bot, speed):
    if random.random() > 0.5:
        ran_theta = np.pi / 2
    else:
        ran_theta = -np.pi / 2

    sensor_mat = [[1, 0, 0, ran_theta],
                  [1, 1, 0, -np.pi / 2],
                  [1, 0, 1, np.pi / 2],
                  [1, 1, 1, np.pi]]

    for i in range(4):
        if front_sensor == sensor_mat[i][0] and left_sensor == sensor_mat[i][1] and right_sensor == sensor_mat[i][2]:
            add_theta = sensor_mat[i][3]
            break
        else:
            add_theta = 0

    resultant_heading = bot.heading + add_theta

    newix = bot.ix + speed * math.cos(resultant_heading)
    newiy = bot.iy + speed * math.sin(resultant_heading)

    return resultant_heading, newix, newiy

def obstacle_avoidance_velocity(left_sensor, front_sensor, right_sensor, bot):
    if random.random() > 0.5:
        ran_theta = np.pi / 2
    else:
        ran_theta = -np.pi / 2

    sensor_mat = [[1, 0, 0, ran_theta],
                  [1, 1, 0, -np.pi / 2],
                  [1, 0, 1, np.pi / 2],
                  [1, 1, 1, np.pi],
                  [0, 1, 0, -np.pi/2],
                  [0, 0, 1, +np.pi/2]]

    for i in range(6):
        if front_sensor == sensor_mat[i][0] and left_sensor == sensor_mat[i][1] and right_sensor == sensor_mat[i][2]:
            add_theta = sensor_mat[i][3]
            break
        else:
            add_theta = 0

    resultant_heading = bot.heading + add_theta

    bot.velocity[0] = 0.05 * math.cos(resultant_heading)
    bot.velocity[1] = 0.05 * math.sin(resultant_heading)

    return bot.velocity

def levy_movement():
    beta = 1
    velocity = [0, 0]
    num = math.gamma(1 + beta) * np.sin(np.pi * beta / 2)
    den = math.gamma((1 + beta) / 2) * beta * 2 ** ((beta - 1) / 2)
    sigma_u = (num / den) ** (1 / beta)
    k = 0  # random.randint(-50, 50)
    u = np.random.normal(k, sigma_u ** 2, 1)
    v = np.random.normal(k, 1, 1)
    step_size = abs(u[0] / (abs(v[0]) ** (1 / beta)))
    randnum = [-1, 1]
    theta = randnum[random.randint(0, 1)] * np.pi * random.random()

    velocity[0]  = step_size * math.cos(theta)
    velocity[1]  = step_size * math.sin(theta)
    # if step_size > 20:
    #     step_size = 20
    # if step_size < 1.2:
    #     step_size = 1.2

    return velocity, step_size


def global_variable(distance, j, rpos, grp):
    costheta = 0
    sintheta = 0
    total_dist = 0
    velocity = [0, 0]
    for k in [h for h in grp if h != j]:
        # print(k)
        theta = math.atan2(rpos[j][1] - rpos[k][1], rpos[j][0] - rpos[k][0])
        costheta = (1 / distance[j][k]) * (math.cos(theta)) + costheta
        sintheta = (1 / distance[j][k]) * (math.sin(theta)) + sintheta
        total_dist = (1 / distance[j][k]) + total_dist

    costheta = ((costheta / total_dist))
    sintheta = ((sintheta / total_dist))
    restheta = math.atan2(sintheta, costheta)
    velocity[0] = total_dist * math.cos(restheta)
    velocity[1] = total_dist * math.sin(restheta)

    return velocity

class motionplanning_rogue:

    def __init__(self, j, rx, ry):
        self.levy_heading = 0
        self.step_size = 0
        self.state = 3
        self.check = 0
        self.current_distance = 0
        self.ix, self.iy = rx[j], ry[j]
        self.heading = 0.00
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

        self.left_sensor, self.front_sensor, self.right_sensor = sensor_reading(current_grid, self,
                                                                                                   pixel,
                                                                                                   self.sensor_range)

        if self.left_sensor + self.front_sensor + self.right_sensor > 0:
            self.state = 1
        else:
            self.state = 3

        if self.state == 1:
            self.velocity = obstacle_avoidance_velocity(self.left_sensor, self.front_sensor, self.right_sensor, self)
            self.step_size = 0
        else:
            if (-self.current_distance + self.step_size) < 0.05:
                self.lbest_velocity, self.step_size = levy_movement()
                self.current_distance = 0


            self.velocity = self.lbest_velocity

        self.heading = math.atan2(self.velocity[1], self.velocity[0])

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

    def levystepsize(self):
        if abs(self.current_distance - self.step_size) < 0.2:
            self.levy_heading, self.step_size = levyheading()
            # print(step_size[j])
            self.current_distance = 0

        return self.levy_heading, self.step_size, self.current_distance
