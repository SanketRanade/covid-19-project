#include<stdio.h>
#include<stdlib.h>
#include<time.h>
 
#define MAX_VERTICES 10000
#define MAX_EDGES 3000

int size = 0;						//for storing size of binary heap
int status[MAX_VERTICES];			//0 means susceptible, 1 means infected, 2 means recovered
int pred_inf_time[MAX_VERTICES];	//for storing time instance at which node gets infected
int recovery_time[MAX_VERTICES]; 	//for storing time instance at which node gets reflected

int find_min(int a, int b, int c)		//function for finding minimum of 3 numbers
{
	if (a <= b && a <= c)
	{
		return a;
	}
	if (b < a && b < c)
	{
		return b;
	}
	return c;
}

typedef struct node_queue		//struct for implementation of priority queue using binary heap
{
	int action;					//for storing action (transmit = 1, recover = 2)
	int time;					//for storing time at which this process happens
	int node_number;			//number in the graph
}node_queue;

int parent(int i)			//for finding parent of node in binary heap
{
	return (i-1)/2;
}

int leftchild(int i)		//for finding left child of node in binary heap
{
	return (2*i + 1);
}

int rightchild(int i)		//for finding right child of node in binary heap
{
	return (2*i + 2);
}

void swap(node_queue *x, node_queue *y)			//function for swapping 2 nodes in binary heap
{ 
	node_queue temp = *x;
	*x = *y;
	*y = temp;
}

void shiftup(int i, node_queue q[])				//function for moving particular element up in binary heap to satisfy min-heap property
{
	while (i != 0 && q[parent(i)].time > q[i].time)		//stop when your reach at root or when min-heap property is satisfied
	{
		swap(&q[i], &q[parent(i)]);				//move element up in heap
		i = parent(i);						//modify index i
	}
}

void shiftdown(int i, node_queue q[])			//function for moving particular element down in binary heap to satisfy min-heap property
{
	int l = leftchild(i);				//left child 
	int r = rightchild(i);				//right child
	int smallest = i;
	if (l < size && q[l].time < q[smallest].time)
	{
		smallest = l;
	}
	if (r < size && q[r].time < q[smallest].time)
	{
		smallest = r;
	}
	if (smallest != i)
	{
		swap(&q[i], &q[smallest]);			//swap with maximum of left and right child if min-heap property is violated
		shiftdown(smallest, q);				//call recursively in lower subtree
	}
}

void insert(node_queue q[], node_queue k)		//function for inserting node in heap
{
	if (size == MAX_VERTICES)
	{
		printf("Overfflow: Could not insert key\n");
		return;
	}
	int i = size;				//first insert new node in last position
	q[i] = k;					
	size++;						//update size
	shiftup(i, q);				//fix min-heap property if violated by calling shiftup function
}

node_queue extract_min(node_queue q[])
{
	if (size <= 0)
	{
		node_queue k = {-1, -1, -1};
        printf("Error, no elements are present\n");
		return k;
	}
	if (size == 1)			//if only one element is present in heap
	{
		size--;
		return q[0];
	}
	node_queue root = q[0];			//store root node in variable
	q[0] = q[size-1];				//put last element in heap at top position in heap
	size--;							//update size
	shiftdown(0, q);				//fix min-heap property if violated by calling shiftdown function
	return root;					//return original root node stored in variable
}

//helper function for finding index of node with particular node number
//This function will be used for changing priority
int is_in_queue(node_queue q[], int node_no, int new_time)			
{
	int ans = -1;
	for (int i = 0; i < size; i++)
	{
		if (node_no == q[i].node_number)
		{
			if (new_time < q[i].time)
			{
				return i;			//return only if above 2 conditions are satisfied
			}
		}
	}
	return ans;			//if element is not found, return -1
}

void change_priority(node_queue q[], int i, int new_time)	//function for changing priority of particular node in binary heap
{
	q[i].time = new_time;		//assign new value
	shiftup(i, q);				//as new value has higher priority, call shiftup function to restore min-heap property if violated
}

typedef struct AdjListNode
{
	int data;						//stores value of connected nodes(neighbours) in graph
	struct AdjListNode* next;		//pointer for pointing to next node
}AdjListNode;

typedef struct AdjList
{
	AdjListNode* head;			//stores head of adjacency list for each node
}AdjList;

typedef struct graph
{
	int N;					//stores size of total nodes in graph
	AdjList* arr;			//pointer to array of adjacency list head pointers
}graph;

AdjListNode* NewNode(int x)					//utility function for adding new node as neighbour
{
	AdjListNode* temp = (AdjListNode *)malloc(sizeof(AdjListNode));			//allocate memory
	temp->data = x;					//assign value
	temp->next = NULL;
	return temp;
}

graph* create_graph(int x)				//function for creating an empty graph
{
	graph* new_graph = (graph *)malloc(sizeof(graph));				//allocate memory of appropriate size
	new_graph->N = x;												//assign value
	new_graph->arr = (AdjList *)malloc(x * sizeof(AdjList));		//allocate memory for array of adjacency list head pointers
	for (int i = 0; i < x; i++)
	{
		((new_graph->arr)+i)->head = NULL;			//assign values to NULL
	}
	return new_graph;
}

graph* AddEdge(graph* gr, int source, int dest)			//function for adding new edge in graph(for assigning neighbours)
{
	AdjListNode* temp = NewNode(dest);				//call NewNode function
	temp->next = ((gr->arr)+source)->head;			//insert in front of adjacency list of source node
	((gr->arr)+source)->head = temp;
	
	//as this is undirected graph both ways links are necessary
	temp = NewNode(source);							//call NewNode function
	temp->next = ((gr->arr)+dest)->head;			//insert in front of adjacency list of dest node
	((gr->arr)+dest)->head = temp;
	return gr;
}

void PrintGraph(graph* gr)					//function for printing graph
{
	for (int i = 0; i < gr->N; i++)
	{
		AdjListNode* ptr = ((gr->arr)+i)->head;
		printf("Adjacency list of vertex %d\nhead", i);
		while (ptr != NULL)
		{
			printf("-> %d ",ptr->data);
			ptr = ptr->next;
		}
		printf("\n");	
	}	
}

void fast_SIR(graph* gr, double taw, double gamma, int initial_infecteds[], int size_initial_infecteds, int tmax);
void process_trans_SIR(graph* gr, node_queue event, node_queue q[], double taw, double gamma, int t, int tmax);
void find_trans_SIR(node_queue q[], double taw, int t, int source_rec_time, AdjListNode* target, int tmax);
void process_rec_SIR(node_queue node, int t);
void display_SIR_numbers();

int coin_prob(int t, double prob)			//function for tossing a biased coin with probability prob
{
	int i = 0;
	int prob_convrt = (int)(prob*100);
	int number = 100;
	while (number > prob_convrt)			//generate random number until you get required result
	{
		number = rand() % 100 + 1;
		i++;
	}
	int retrn_time = t + i;			//time which we are going to return
	return retrn_time;
}

void process_trans_SIR(graph* gr, node_queue event, node_queue q[], double taw, double gamma, int t, int tmax)
{
	int num1 = event.node_number;
	status[num1] = 1;						//update status of that particular node (status = 1 means node has infected)
	int rec_time = coin_prob(t, gamma);		//calculate day at which this node will recover
	if (rec_time < tmax)					//if recovery time is less than maximum time
	{
        node_queue k = {2, rec_time, num1};		//store required information in new node (action = 2 means recover)
		insert(q, k);								//insert this node in priority queue(binary heap)
	}
	AdjListNode* ptr = ((gr->arr)+num1)->head;		//go to the first node in neighbours of this node in graph
	while (ptr != NULL)								//go through all neighbours
	{
		find_trans_SIR(q, taw, t, rec_time, ptr, tmax);			//call find_trans_SIR function for each neoghbours
		ptr = ptr->next;									//go to the next neighbour
	}
}

void find_trans_SIR(node_queue q[], double taw, int t, int source_rec_time, AdjListNode* target, int tmax)
{
	int num2 = target->data;
	if (status[num2] == 0)					//if target is node is susceptible
	{
		int inf_time = coin_prob(t, taw);					//calculate day at which this node will get infected
		if (inf_time < find_min(source_rec_time, pred_inf_time[num2], tmax))	//node can get infected from this source only if this condition is satisfited
		{
			//this means this node was already getting infected by other source, but infection day from that node was late
			if (pred_inf_time[num2] != __INT_MAX__)	
			{
				int j = is_in_queue(q, num2, inf_time);			//find index of target node in priority queue
				change_priority(q, j, inf_time);				//update priority of target node in priority queue
			}
			//this means this node hasn't infected from any other source 
			else
			{
				node_queue k = {1, inf_time, num2};		//store required information in new node (action = 1 means transmit)
				insert(q, k);    							//insert this node in priority queue(binary heap)
			}
			pred_inf_time[num2] = inf_time;			//update infection time of this node as inf_time
		}
	}
}

void process_rec_SIR(node_queue node, int t)
{
	int num = node.node_number;
	status[num] = 2;					//update status of particular node (status = 2 means node has recovered)
	recovery_time[num] = t;				//update recovery time of that node in recovery_time array
}

void fast_SIR(graph* gr, double taw, double gamma, int initial_infecteds[], int size_initial_infecteds, int tmax)
{
	node_queue q[MAX_VERTICES];						//assigning variable for binary min-heap
	for (int i = 0; i < MAX_VERTICES; i++)
	{
		status[i] = 0;							//assign status to all nodes (0 means susceptible)
		pred_inf_time[i] = __INT_MAX__;			//this is predefined infection time
		recovery_time[i] = -1;					//this is time at which a particular node will recover (-1 means person is not infected)
	}
	
	for (int i = 0; i < size_initial_infecteds; i++)				//loop through initial infected nodes
	{
		node_queue k = {1, 0, initial_infecteds[i]};			//store required information in new node (action = 1 means transmit)
		insert(q, k);										//insert this node in priority queue(binary heap)
		pred_inf_time[initial_infecteds[i]] = 0;			//update infection time of this node as 0
	}

	while (size > 0)						//loop until priority queue becomes empty
	{
		node_queue event = extract_min(q);			//remove element with highest priority
		int num = event.node_number;
		int current_time = event.time;
		if (event.action == 1)					//if action is transmit
		{
			if (status[num] == 0)				//check whether node is susceptible or not
			{
				process_trans_SIR(gr, event, q, taw, gamma, current_time, tmax);		//call process_trand_SIR function
			}
		}
		else if (event.action == 2)			//if action is recover
		{
			process_rec_SIR(event, current_time);		//call process_rec_SIR function
		}

		if (size > 0)
		{
			if (q[0].time != current_time)			//print after doing all actions for particular time
			{
				int temp_time = current_time;
				while (temp_time != q[0].time)
				{
					printf("day %d: ", temp_time);		//prints day number
					display_SIR_numbers();				//print numbers of each category
					printf("\n");
					temp_time++;
				}
			}			
		}
		if (size == 0)					//printing after all actions are done
		{
			printf("day %d: ", current_time);
			display_SIR_numbers();					//print numbers of each category
			printf("\n");
		}	
	}
}

void display_SIR_numbers()				//function for displaying numbers of each category
{
	int cnt0 = 0;
	int cnt1 = 0;
	int cnt2 = 0;
	for (int i = 0; i < MAX_VERTICES; i++)
	{
		if (status[i] == 2)				//status = 2 means node is recovered
		{
			cnt2++;
		}
		else if (status[i] == 1)		//status = 1 means node is infected
		{
			cnt1++;
		}
		else							//status = 1 means node is susceptible
		{
			cnt0++;
		}
	}
	printf("total recovered is %d, total infected is %d, total non-infected is %d\n", cnt2, cnt1, cnt0);
}

int main()
{
	srand(time(NULL));
	int NumVertices = MAX_VERTICES;	
	graph* grp = create_graph(NumVertices);					//create an empty graph

	int NumEdges = 1 + rand() % MAX_EDGES;					//randomly generate number of edges in graph
	while ((NumVertices)*(NumVertices-1)/2 < NumEdges)		//check whether these much edges are possible or not in graph
	{
		NumEdges = 1 + rand() % MAX_EDGES;
	}
	
	int edge_counter = 0;
	for (edge_counter = 0; edge_counter < NumEdges; edge_counter++)
	{
		int flag = 0;
		int vertex1 = rand() % NumVertices;				//randomly generate 2 vertices
		int vertex2 = rand() % NumVertices;
		if (vertex1 == vertex2)					//if both vertices generated are same, this edge is not possible
		{
			flag = 1;
			edge_counter--;
		}
		else				//check whether this edge already exists or not
		{
			for(AdjListNode* ptr = ((grp->arr)+vertex1)->head; ptr != NULL; ptr = ptr->next)	//loop through all neighbours of node
			{
				if (ptr->data == vertex2)			//if already edge is created (if already a neighbour)
				{
					flag = 1;
					edge_counter--;
					break;
				}
			}
		}
		if (flag == 0)
		{
			grp = AddEdge(grp, vertex1, vertex2);		//add edge between 2 vertices
		}
	}

	//PrintGraph(grp);

	printf("total vertices in graph = %d, total edges = %d\n\n",NumVertices, NumEdges);

	int size_initial_infecteds = 1 + rand() % MAX_VERTICES;			//randomly generate size of initially infected people

	int initial_infecteds_list[size_initial_infecteds];			//create array of required size
	int k = 0;
	for (k = 0; k < size_initial_infecteds; k++)
	{
		int fl = 0;
		int r = rand() % MAX_VERTICES;			//randolmy mark node as initially infected
		for (int t = 0; t < k; t++)				//check whether this node is already present in array or not
		{
			if (r == initial_infecteds_list[t])
			{
				fl = 1;
				k--;
				break;
			}
		}
		if (fl == 0)			//if node is not present in array of initially infected people
		{
			initial_infecteds_list[k] = r;		//put that node number in array
		}
	}

	//call fast_SIR function with taw = 0.5, gamma = 0.2, maximum time = 300 days
	fast_SIR(grp, 0.5, 0.2, initial_infecteds_list, size_initial_infecteds, 300);

	printf("total vertices in graph = %d, total edges = %d\n",NumVertices, NumEdges);
	printf("initial infecteds = %d\n", size_initial_infecteds);

	return 0;
}
