/*
*
Haibo Ni - haibo.ni02@gmail.com

MPI toolbox for tissue cell-to-cell connectivity and MPI communication
local to global message passing.
* code modified based on https://github.com/AlexeyMalkhanov/Cardiac_demo/blob/master/heart_demo.cpp, with significant modifications and enhancements
* copyrights of the original code (and thus this file) see https://github.com/AlexeyMalkhanov/Cardiac_demo
*/



#ifndef TOPOLOGY_H
#define TOPOLOGY_H
#include <mpi.h>
#include <map>
#include <vector>
#include <iostream>
// #include <iostream>
#include <algorithm>
// #include <vector>
#include <iterator>
// #define DEBUG_TOP (1)
class Topology
{
public:
	Topology(MPI_Comm comm, int cell_num) {
		// int rank, size;
		MPI_Comm_rank(comm, &rank);
		MPI_Comm_size(comm, &process_num);
		char processor_name[MPI_MAX_PROCESSOR_NAME];
		mpi_comm = comm;

		total_cell_num = cell_num;
		local_cell_num = total_cell_num / process_num;

		// if (total_cell_num % process_num)
		if (rank < total_cell_num % process_num) {
			local_cell_num = local_cell_num + 1;
		}
		// printf("Hybrid: Hello from thread from process %d out of %d on %s, cell num is %d\n",
		//        rank, process_num, processor_name, local_cell_num);


		local_cell_num_vec  = new int [process_num];
		offset_cell_num_vec = new int  [process_num];

		for (int i = 0; i < process_num; ++i)
		{
			local_cell_num_vec[i] = total_cell_num / process_num;
			if (i < total_cell_num % process_num) {
				local_cell_num_vec[i] ++ ;
			}
		}

		starting_index = 0;

		offset_cell_num_vec[0] = 0;


		// assign cells to each processor/rank cpus
		for (int i = 0; i < process_num; ++i)
		{
			if (i < rank) {
				starting_index += local_cell_num_vec[i];

			}
		}

		for (int i = 1; i < process_num; ++i)
		{
			offset_cell_num_vec[i] = offset_cell_num_vec[i - 1] + local_cell_num_vec[i - 1];
		}

		if (rank == 2)
			for (int i = 0; i < process_num; ++i)
			{
				printf("***********cell %d, offset %d \n",
				       local_cell_num_vec[i], offset_cell_num_vec[i]);
			}


		for (int jj = 0; jj < process_num; ++jj)
		{


			std::vector<int> v;

			for (int i = offset_cell_num_vec[jj]; i < offset_cell_num_vec[jj] + local_cell_num_vec[jj]; ++i)
			{
				v.push_back(i);
			}
			total_rank_id_map.insert(std::make_pair(jj, v));
		}


		// printf("Hybrid: Hello from thread from process %d out of %d on %s, cell num is %d, starting at %d\n",
		// rank, process_num, processor_name, local_cell_num, starting_index);

	};


	void setup_dependecies(int * nbd, const int ID_nbd) {

		// looking for cells with id not listed in the current rank processor
		// these cells will be the gohst cells;
		for (int i = 0; i < local_cell_num; ++i)
		{

			for (int ii = 1; ii < ID_nbd; ++ii)
			{
				int ID = nbd[((i + starting_index) * ID_nbd) + ii];

				if (ID < starting_index or ID >= starting_index + local_cell_num)
				{
					non_local_ids.push_back(ID);
				}
			}
		}


		// std::cout << " ****************** rank non_local_ids" << rank << " " ;
		// for (auto &c : non_local_ids)
		// {
		// 	std::cout << c << " ";
		// }
		// std::cout << std::endl;
		// now we need to know where are the ghost cells (in which rank are the cells processed)
		for (int jj = 0; jj < process_num; ++jj)
		{

			std::vector<int> v;

			for (int i = offset_cell_num_vec[jj]; i < offset_cell_num_vec[jj] + local_cell_num_vec[jj]; ++i)
			{

				for (int ii = 1; ii < ID_nbd; ++ii)
				{
					int ID = nbd[((i/* + offset_cell_num_vec[jj]*/) * ID_nbd) + ii];

					if (ID < offset_cell_num_vec[jj] or ID >= offset_cell_num_vec[jj] + local_cell_num_vec[jj])
					{
						v.push_back(ID);
						// all_non_local_ids_map[jj].push_back(ID);
					}
				}

			}
			all_non_local_ids_map.insert(std::make_pair(jj, v));
			v.clear();
		}


#ifdef DEBUG_TOP
		if (rank == 1) {
			for (const auto &pair : all_non_local_ids_map) {
				std::cout << " **$$$$$$$$$$$$$$ first " << pair.first << " " ;
				for (auto &c : pair.second)
				{
					std::cout << c << " ";
				}
				std::cout << std::endl;

			}

		}
#endif


		// create non local dependency map:

		for (auto &c : non_local_ids)
		{

			for (const auto &pair : total_rank_id_map) {
				// std::vector<int> v;
				if (std::find(std::begin(pair.second), std::end(pair.second), c)  != std::end(pair.second))   // -> c in pair.second, if find c in the vec, ghost cell in the pair.first rank nodes
					// v.push_back()
					non_loc_dependecies[pair.first].push_back(c);

			}
			// if (v.size() > 0)
		}

#ifdef DEBUG_TOP

		if (rank == 1) {
			for (const auto &pair : non_loc_dependecies) {
				std::cout << " &&&&&&&&&&&&&&&& first " << pair.first << " " ;
				for (auto &c : pair.second)
				{
					std::cout << c << " ";
				}
				std::cout << std::endl;

			}
		}
#endif

		// create non_local_dependee map: loop through all other ranks, find all ghost cells that originates from the current rank node block
		for (auto &pair : all_non_local_ids_map)
		{

			// total_rank_id_map[rank] list all cell id computed in the current mpi rank

			for (const auto & id : pair.second)  {
				if (std::find(std::begin(total_rank_id_map[rank]), std::end(total_rank_id_map[rank]), id)  != std::end(total_rank_id_map[rank]))
					non_loc_dependees[pair.first].push_back(id);
			}
		}


#ifdef DEBUG_TOP
		if (rank == 2) {
			for (const auto &pair : non_loc_dependees) {
				std::cout << " MMMMMMMMMMM first " << pair.first << " " ;
				for (auto &c : pair.second)
				{
					std::cout << c << " ";
				}
				std::cout << std::endl;

			}
		}
#endif

		// set up non_loc_dependees bufs
		for (auto &n : non_loc_dependees) {
			std::vector<double> buf(n.second.size());

			for (int i = 0; i < buf.size(); ++i)
			{
				buf[i] = (double) n.second[i] + 0.1 * (i + 1);
				std::cout << buf[i] << " ";
			}
			non_loc_dependees_bufs[n.first] = buf;
		}

		// set up non_loc_dependees bufs
		for (auto &n : non_loc_dependecies) {
			std::vector<double> buf(n.second.size());
			non_loc_dependecies_bufs[n.first] = buf;
		}


		sreqs.resize(non_loc_dependees.size());
		rreqs.resize(non_loc_dependecies.size());
		sstats.resize(non_loc_dependees.size());
		rstats.resize(non_loc_dependecies.size());
	}



	inline void start_sends() {
		int count = 0;
		for (auto &d : non_loc_dependees ) {
			int remote_rank = d.first;
			double *buf = non_loc_dependees_bufs[remote_rank].data();
			int size = non_loc_dependees_bufs[remote_rank].size();
			MPI_Isend(buf, size, MPI_DOUBLE, remote_rank, 123,
			          mpi_comm, &sreqs[count++]);

		}
	}

	inline void start_recvs() {
		int count = 0;
		int offset = 0;
		char *recv_buf = (char*)non_local_values.data();
		for (auto &d : non_loc_dependecies ) {
			int remote_rank = d.first;
			int size = d.second.size() * sizeof(double);
			double *buf = non_loc_dependecies_bufs[remote_rank].data();
			MPI_Irecv(buf, size + 10, MPI_DOUBLE,
			          remote_rank, 123, mpi_comm,
			          &rreqs[count++]);
			// offset += size;
		}
	}

	inline void wait_sends() {
		MPI_Waitall(sreqs.size(), sreqs.data(), sstats.data());
	}
	inline void wait_recvs() {
		MPI_Waitall(rreqs.size(), rreqs.data(), rstats.data());
	}

	inline void exchange() {

		MPI_Barrier(mpi_comm);
		start_recvs();
		start_sends();
		wait_sends();
		wait_recvs();

#ifdef DEBUG_TOP

		if (rank == 1) {
			for (const auto &pair : non_loc_dependecies_bufs) {
				std::cout << " non_loc_dependecies_bufs first " << pair.first << " " ;
				for (auto &c : pair.second)
				{
					std::cout << c << " ";
				}
				std::cout << std::endl;

			}
		}
#endif

	}


	inline void set_up_send_bufs(double * v_new) {
		// SINGLETHREAD
		for (auto &d : non_loc_dependees) {
			int remote_rank = d.first;
			std::vector<int> &ids_to_pack = d.second;
			std::vector<double> &buf = non_loc_dependees_bufs[remote_rank];
			for (int i = 0; i < ids_to_pack.size(); i++) {
				int j = ids_to_pack[i];
				buf[i] = v_new[j - starting_index]; // v_new is local to processor and j is the global index;
				// buf[i] = cnodes[j].state[DynamicalSystem::COUPLING_VAR_ID];
				// if (rk4_id == 3)
				// 	buf[i] += cnodes[j].rk4[2][DynamicalSystem::COUPLING_VAR_ID] / 2.0;
				// else if (rk4_id > 0)
				// 	buf[i] += cnodes[j].rk4[rk4_id - 1][DynamicalSystem::COUPLING_VAR_ID];
			}
		}
	}


	inline void map_recvs_bufs(double * v_global) {
		// SINGLETHREAD
		for (auto &d : non_loc_dependecies) {
			int remote_rank = d.first;
			std::vector<int> &ids_to_pack = d.second;
			std::vector<double> &buf = non_loc_dependecies_bufs[remote_rank];
			for (int i = 0; i < ids_to_pack.size(); i++) {
				int j = ids_to_pack[i];
				v_global[j] = buf[i];
				// buf[i] = v_new[j - starting_index]; // v_new is local to p
			}
		}
	}

	inline void map_local_and_ghost_to_global(double * v_new, double * v_global) {
		// SINGLETHREAD
		memcpy ( &v_global[starting_index], v_new, local_cell_num * sizeof(double) ); // fast copy using memcpy
	}


	inline void set_ghost_cells(double * v_new, double * v_global) {

		if (process_num == 1) {

			// std::cout << "sasssss" << std::endl;
			// if only one processor used, then just copy everything here.
			memcpy ( v_global, v_new, local_cell_num * sizeof(double));
		} else {

			// int error = MPI_Gatherv(v_new, local_cell_num, MPI_DOUBLE, v_global, local_cell_num_vec, offset_cell_num_vec, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			// error = MPI_Bcast(v_global, total_cell_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			// single threaded
			set_up_send_bufs(v_new);
			exchange();

			map_recvs_bufs(v_global);

			map_local_and_ghost_to_global(v_new, v_global);
		}
	}

	~Topology() {

		if (local_cell_num_vec) {
			delete [] local_cell_num_vec;
		}
		if (offset_cell_num_vec) {
			delete [] offset_cell_num_vec;
		}
	};

	int rank;
	int process_num;
	int total_cell_num;
	int local_cell_num;
	int starting_index;

	int* local_cell_num_vec;
	int* offset_cell_num_vec;

	MPI_Comm mpi_comm;

	std::vector<double> non_local_values;
	std::vector<int> non_local_ids;
	std::vector<int> local_ids;
	std::map<int, std::vector<int> > total_rank_id_map;
	std::map<int, std::vector<int> > all_non_local_ids_map;
	std::map<int, std::vector<int> > non_loc_dependecies;
	std::map<int, std::vector<int> > non_loc_dependees;
	std::map<int, std::vector<double> > non_loc_dependees_bufs;
	std::map<int, std::vector<double> > non_loc_dependecies_bufs;
	std::vector<MPI_Request> sreqs;
	std::vector<MPI_Request> rreqs;
	std::vector<MPI_Status> sstats;
	std::vector<MPI_Status> rstats;

};


#endif