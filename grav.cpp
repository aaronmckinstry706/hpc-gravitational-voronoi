#include <stdio.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string.h>
#include <string>
#include <iostream>
#include <sstream>
#include <ctime>
#include <unistd.h>
#include <pthread.h>

#define GRID_SIZE (1000)
#define NUM_THREADS (2)
#define NO_THREAD_AVAILABLE (-1)
#define CIRCLE_INCREMENT (10)

using namespace std;

struct grav_point {
	int x;
	int y;
};


struct move {
	int x;
	int y;
	int player;
	int score;
	
	std::string move_string();	
};

struct thread_params {
    struct move* m_ptr;//ptr to location at which we will store return value
    struct move a;
    int player;
    bool* inUse_ptr;
    
    thread_params() {}
    
    thread_params(struct move* m, struct move aa, int p)
        : m_ptr(m), a(aa), player(p) {}
};

class thread_pool {
    pthread_t threads[NUM_THREADS];
    bool inUse[NUM_THREADS];
    thread_params thread_params_list[NUM_THREADS];
    
    int getNextAvailableThread() {
        for (int i = 0; i < NUM_THREADS; ++i) {
            
            if (!inUse[i])
                return i;
        }
        return -1;
    }
    
public:
    
    thread_pool() {
        for (int i = 0; i < NUM_THREADS; ++i)
            inUse[i] = false;
    }
    
    void call(void* (*f)(void*), thread_params tps, int waitTimeMicroseconds) {
        int availableThread = getNextAvailableThread();
        while (availableThread == NO_THREAD_AVAILABLE) {
            
            usleep(waitTimeMicroseconds);
            availableThread = getNextAvailableThread();
        }
        thread_params_list[availableThread] = tps;
        thread_params_list[availableThread].inUse_ptr = &inUse[availableThread];
        inUse[availableThread] = true;
        int rc = pthread_create(&threads[availableThread], NULL, 
                            f, (void*)&thread_params_list[availableThread]);
    }
    
    void joinAll() {
        
        for (int i = 0; i < NUM_THREADS; ++i) {
            if (inUse[i]) {
                pthread_join(threads[i],NULL);
            }
        }
    }
};


int distances[GRID_SIZE][GRID_SIZE];
double pull[2*GRID_SIZE*GRID_SIZE];

double grid[GRID_SIZE][GRID_SIZE];
int stones[GRID_SIZE][GRID_SIZE];

// points 66 closest to stone
vector<struct grav_point> grav_points;


int player = -1;
int NUMBER_STONES = 0;

std::string my_itoa(int number) {
	stringstream oss;
	oss << number;
	return oss.str();
}

std::string move::move_string() {
	std::string s = my_itoa(this->x) + " " +my_itoa(this->y);
	return s;
}

double get_distance(int x1, int y1, int x2, int y2) {
	int sum = abs(x1-x2) << 1;
	sum += abs(y1-y2) << 1;
	return sqrt(sum);
}

void build_distance_table() {

	for (int i = 0; i < GRID_SIZE; i++) {
		for (int j = 0; j < GRID_SIZE; j++) {
			distances[i][j] = i*i + j*j;
		}
	}
}

void build_pull_table() {

	for (int i = 0; i < 2*GRID_SIZE*GRID_SIZE; i++) {
		pull[i] = 1.0/i;
	}
}

void initialize_grid() {
	memset(&grid, 0, GRID_SIZE*GRID_SIZE*sizeof(double));
}


// TODO: make actual circle instead of psuedocircle
void generate_points() {

	struct grav_point gp;

	// going CCW, starting on right, going up
	int current_x = 0;
	int current_y = 0;
	for (int i = 0; i <= 66; i+=CIRCLE_INCREMENT) {
		
		while (sqrt(distances[abs(current_x)][abs(current_y)]) < 66) {
			current_x++;
		}
		gp.x = current_x;
		gp.y = current_y;
		current_y++;
		current_x = 0;
		grav_points.push_back(gp);
	}

	// starting on top, going left
	current_x = 0;
	current_y = 0;
	for (int i = 0; i <= 66; i+=CIRCLE_INCREMENT) {
		while (sqrt(distances[abs(current_x)][abs(current_y)]) < 66) {
			current_y++;
		}
		gp.x = current_x;
		gp.y = current_y;
		current_y = 0;
		current_x--;
		grav_points.push_back(gp);
	}
	
	// starting on left, going down
	current_x = 0;
	current_y = 0;
	for (int i = 0; i <= 66; i+=CIRCLE_INCREMENT) {
		while (sqrt(distances[abs(current_x)][abs(current_y)]) < 66) {
			current_x--;
		}
		gp.x = current_x;
		gp.y = current_y;
		current_y--;
		current_x = 0;
		grav_points.push_back(gp);
	}
	
	// starting on left, going down
	current_x = 0;
	current_y = 0;
	for (int i = 0; i <= 66; i+=CIRCLE_INCREMENT) {
		while (sqrt(distances[abs(current_x)][abs(current_y)]) < 66) {
			current_y--;
		}
		gp.x = current_x;
		gp.y = current_y;
		current_y = 0;
		current_x++;
		grav_points.push_back(gp);
	}
}

// apply last move to stones array, update grid with last move
void update_grid(vector<struct move> moves) {
	
	if (moves.size() == 0)
		return;
	struct move m = moves.back();
	
	stones[m.x][m.y] = m.player;
//	printf("at %d %d -> %d\n", m.x, m.y, m.player);
	
	if (m.player==player) {
		for (int i = 0; i < GRID_SIZE; i++) {
			for (int j = 0; j < GRID_SIZE; j++) {
				grid[i][j] += pull[distances[abs(m.x - i)][abs(m.y - j)]];
			}
		}
	} else{
		for (int i = 0; i < GRID_SIZE; i++) {
			for (int j = 0; j < GRID_SIZE; j++) {
				grid[i][j] -= pull[distances[abs(m.x - i)][abs(m.y - j)]];
			}
		}
	}
	
	for (int i = 0; i < moves.size(); i++) {
		if (moves[i].player == player)
			grid[moves[i].x][moves[i].y] = 1;
		else 
			grid[moves[i].x][moves[i].y] = -1;
	}
}

// generate score based on existing board and new move m (score is # points owned)
int score_new_move(struct move m) {

	// double new_grid[GRID_SIZE][GRID_SIZE];
	// memcpy(&new_grid, &grid, GRID_SIZE*GRID_SIZE*sizeof(double));

	int score = 0;
	for (int i = 0; i < GRID_SIZE; i++) {
		for (int j = 0; j < GRID_SIZE; j++) {
			// new_grid[i][j] += pull[distances[abs(m.x - i)][abs(m.y - j)]];
			if (grid[i][j] + pull[distances[abs(m.x - i)][abs(m.y - j)]] > 0.0)
				score++;
		}
	}
	return score;
}

// check if theres a stone within 66 of suggested move m
bool check_nearby(struct move m) {

	for (int i = max(0, m.x-66); i <= min(999, m.x+66); i++) {
		for (int j = max(0, m.y-66); j <= min(999, m.y+66); j++) {
			if (stones[i][j] != 0) {
				if (sqrt(distances[abs(m.x - i)][abs(m.y - j)]) < 66)
					return false;
			}
		}
	}
	return true;
}

// check if move is out of bounds or has nearby stone
bool is_valid(struct move m) {

	if (m.x < 0 || m.x >= GRID_SIZE || m.y < 0 || m.y >= GRID_SIZE)
		return false;
	return check_nearby(m);
}

// get points around stone
vector<struct grav_point> get_points(struct move a) {

	vector<struct grav_point> gps;
	struct grav_point tmp;
	for (int i = 0; i < grav_points.size(); i++) {
		tmp.x = grav_points[i].x + a.x;
		tmp.y = grav_points[i].y + a.y;
		gps.push_back(tmp);
	}
	return gps;
}

// for stone a, return best valid move around it
void *find_best_around_stone(void *params)
{
    
    thread_params* threadParams = (thread_params*)params;
    struct move a = threadParams->a;
    int player = threadParams->player;
    
    
	int highest = 0;
	struct move m;
	m.x = 0;
	m.y = 0;
	struct move current_move;
	int score;
	vector<struct grav_point> gps = get_points(a);
	
//	printf("gps is %d\n", gps.size());
	
	for (int i = 0; i < gps.size(); i++) {
		current_move.x = gps[i].x;
		current_move.y = gps[i].y;
		current_move.player = player;
		if (!is_valid(current_move))
			continue;
		score = score_new_move(current_move);
		current_move.score = score;
		if (score > highest) {
			highest = score;
			m = current_move;
		}
	}
	*(threadParams->m_ptr) = m;
    *(threadParams->inUse_ptr) = false;
    
}



// player is 1 or 2 (finding best move as player <player>)
struct move find_best_around(vector<struct move> moves, int p) {
	int highest = 0;
	struct move m;
	m.x = 0;
	m.y = 0;
	struct move current_move;
	
	struct move best[moves.size()];
	
	clock_t time1 = clock();
	
	// get best per stone
    thread_pool pool;
	for (int i = 0; i < moves.size(); i++) {
		if (moves[i].player != p) {
            
            thread_params tps(&best[i], moves[i], p);
            pool.call(find_best_around_stone, tps, 100000);
		}
	}
    pool.joinAll();
		
	// get best of around all stones
	for (int i = 0; i < moves.size(); i++) {
		if (moves[i].player != p) {
			if (best[i].score > highest) {
				m = best[i];
				highest = best[i].score;
			}
		}
	}
	return m;
}



struct move pick_next_move(vector<struct move> moves) {

	if (moves.size() == 0) {
		struct move m;
		m.x = 500;
		m.y = 500;
		return m;
	}
	
	if (player == 1) {
		struct move m;
		m.x = 0;
		m.y = 0;
		
		if (moves.size() == 2*NUMBER_STONES - 2) {
			m = find_best_around(moves, 1);
		} else {
			struct move last = moves.back();
			m.x = 999-last.x;
			m.y = 999-last.y;
			m.player = 1;
			if (!is_valid(m))  {
				m = find_best_around(moves, 1);
			}
		}
		return m;
	} else {
		return find_best_around(moves, player);
	}
}


int parse_input(char *input, std::vector<struct move> *moves) {
	
	std::string s(input);
	int number;

	std::stringstream ss(s);
	int game_finished = 0;
	ss >> game_finished;
	
	if (game_finished == 1)
		return 1;

	int number_of_moves = 0;
	ss >> number_of_moves;
	for (int i = 0; i < number_of_moves; i++) {
		struct move m;	
		ss >> m.x;
		ss >> m.y;
		ss >> m.player;
		moves->push_back(m);
	}
	return 0;
}

int main(int argc, char *argv[]) {

	int s, port;
	struct sockaddr_in server;

	if (argc < 4) {
		printf("ERROR: call executable with number of stones, players, and port\n");
		return 1;
	}
	NUMBER_STONES = atoi(argv[1]);

	port = atoi(argv[3]);
	
	s = socket(AF_INET, SOCK_STREAM, 0);
	server.sin_addr.s_addr = inet_addr("127.0.0.1");
	server.sin_family = AF_INET;
	server.sin_port = htons(port);
	
	if (connect(s, (struct sockaddr *) &server, sizeof(server)) < 0) {
		printf("error connecting to server\n");
		return 1;
	}

	char buffer[1024];
	while (1) {
		memset(buffer, 0, 1024);
		recv(s, buffer, 1024, 0);
		clock_t time1 = clock();
		std::vector<struct move> moves;
		int status = parse_input(buffer, &moves);
		if (status == 1) {
			break;
		}
		if (moves.size() == 0)
			player = 1;

		if (player == -1) {
			player = moves.size() + 1;
			build_distance_table();
			build_pull_table();
			initialize_grid();
			generate_points();
		}
			
		clock_t time2 = clock();
		update_grid(moves);
		clock_t time3 = clock();
		
		struct move m = pick_next_move(moves);
		clock_t time4 = clock();

		moves.push_back(m);
		update_grid(moves);
		
/*
		// time related stuff
		clock_t time5 = clock();
		printf("parsing move: %f\n", double(time2-time1)/CLOCKS_PER_SEC);
		printf("updating previous move: %f\n", double(time3-time2)/CLOCKS_PER_SEC);
		printf("creating move: %f\n", double(time4-time3)/CLOCKS_PER_SEC);
		printf("updating current move: %f\n", double(time5-time4)/CLOCKS_PER_SEC);
		printf("total move: %f\n", double(time5-time1)/CLOCKS_PER_SEC);
*/
//		printf("%d %d\n", m.x, m.y);
		send(s, m.move_string().c_str(), m.move_string().length(), 0);
	}
	close(s);

	return 0;
}

